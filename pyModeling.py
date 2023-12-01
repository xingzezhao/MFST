#!/usr/bin/env python
# coding: utf-8

import configparser
import MDAnalysis as mda
import numpy as np
import re
import xdrlib
import os
import shutil
import subprocess
import logging
import time
import sys

surfts_type_dict = {}
surfts_ff_dict = {}

def load_surfts_type_mapping():
    surfts_molpath = './database/mol/'
    surftsname2code = {
        'SDS': '000',
        'SDBS' : '001',
        'SLES' : '002',
        'DHSB' : '003',
        'Gemini' : '004',
        'AGO' : '005',
        'AOS' : '006',
        'AES' : '007',
        'AOT' : '008',
        'MES' : '009',
        'OTAB' : '010',
        'DTAB' : '011',
        'CTAB' : '012',
        'CHSB' : '013',
        'MHPB' : '014',
        'BS12' : '015',
        'FAPOE' : '016',
    }
    for i,y in surftsname2code.items():
        surfts_type_dict[i] = surfts_molpath+'s%s.pdb'%y
        surfts_ff_dict[i] = surfts_molpath+'f%s.itp'%y    

def read_config(config_file):
    # 创建配置解析器对象
    config = configparser.ConfigParser()

    # 读取配置文件
    config.read(config_file)

    # 获取配置项的值
    no_surfts = config.getint('pyModeling input', 'no_surfts')
    surfts_type = config.get('pyModeling input', 'surfts_type').split()
    interfc_type = config.get('pyModeling input', 'interfc_type')
    itf_surfts1 = config.get('pyModeling input', 'itf_surfts')
    itf_surfts = [int(value) for value in itf_surfts1.split()]
    waterbox_size = config.get('pyModeling input', 'waterbox_size')
    waterbox_size_list = [float(value) for value in waterbox_size.split()]
    simulationbox_size = config.get('pyModeling input', 'simulationbox_size')
    simulationbox_list = [float(value) for value in simulationbox_size.split()]
    simulationsys_name = config.get('pyModeling input', 'simulationsys_name')
    temp = config.getfloat('pySimulating input', 'Temp')
    packmolpath = config.get('pyModeling input', 'packmolpath')
    WFF = config.get('pySimulating input', 'SimuWatMod')

    return no_surfts, surfts_type, interfc_type,itf_surfts, waterbox_size_list, simulationbox_list, simulationsys_name,temp,packmolpath, WFF

def find_hydrophobic_tail_atoms(pdb_file):
    # 读取PDB文件
    u = mda.Universe(pdb_file)

    # 选择所有碳原子
    carbon_atoms = u.select_atoms('name C*')
    other_atoms = u.select_atoms('not name C* and not name H*')
    # 初始化变量
    max_distance = 0
    farthest_atom_pair = None

    # 遍历所有碳原子对
    for carbon_atom in carbon_atoms:
            for other_atom in other_atoms:
                # 计算两个原子之间的距离
                distance = np.linalg.norm(carbon_atom.position - other_atom.position)

                # 更新最大距离和原子对
                if distance > max_distance:
                    max_distance = distance
                    farthest_atom_pair = (carbon_atom, other_atom)

    if farthest_atom_pair is not None:
        carbon_atom_index = farthest_atom_pair[0].index
        other_atom_index = farthest_atom_pair[1].index
    else:
        print("Error: No hydrophobic tail atoms found.")
        
    return carbon_atom_index, other_atom_index, max_distance

def get_number_of_water(waterbox_size_list,temp):
    dimensions = np.array(waterbox_size_list)
    volume = np.prod(dimensions) * 1e-27
    rho = 1000*(1-0.00021*(temp-277.15))
    m = volume*rho
    n = m*1000/18
    numofwater = round(n*6.022e23)
    return numofwater
    
def generate_packmol_input(no_surfts, surfts_type, interfc_type,itf_surfts, waterbox_size_list, simulationbox_list, 
                           simulationsys_name,hindexlst,tindexlst,distancelst,now):
    simubox = [round(value*10) for value in simulationbox_list]
    center = [round(value*0.5) for value in simubox]
    waterlocd = round((center[2]-waterbox_size_list[2]*5) -5)
    waterlocu = round((center[2]+waterbox_size_list[2]*5) +5)
    waterloc = [0, 0, waterlocd, simubox[0], simubox[1], waterlocu]
    
    

    surflocu = []
    surflocd = []
    anchorlocu = []
    anchorlocd = []
    for dis in distancelst:
        surflocu1 = round((waterlocd) -5)
        surflocd1 = round((surflocu1-dis))
        surflocu.append([0, 0, surflocd1, simubox[0], simubox[1], surflocu1])
        anchorlocu.append([surflocd1+3,surflocu1-3])
        surflocd2 = round((waterlocu) +5)
        surflocu2 = round((surflocd2+dis))
        surflocd.append([0, 0, surflocd2, simubox[0], simubox[1], surflocu2])
        anchorlocd.append([surflocd2+3,surflocu2-3])

    template1 = """
structure {surfs_name}.pdb
    number {num_surfs}
    inside box {surfloc[0]} {surfloc[1]} {surfloc[2]} {surfloc[3]} {surfloc[4]} {surfloc[5]}
atoms {hydrophilic}
    over plane 0 0 1 {hydrophilic_anchor}
end atoms
atoms {hydrophobic}
    below plane 0 0 1 {hydrophobic_anchor}
end atoms
end structure
    """

    template2 = """
structure {surfs_name}.pdb
    number {num_surfs}
    inside box {surfloc[0]} {surfloc[1]} {surfloc[2]} {surfloc[3]} {surfloc[4]} {surfloc[5]}
atoms {hydrophobic}
    over plane 0 0 1 {hydrophobic_anchor}
end atoms
atoms {hydrophilic}
    below plane 0 0 1 {hydrophilic_anchor}
end atoms
end structure
    """

    templateT = """tolerance 1.0
filetype pdb
output {output_filename}.pdb

{stru1}

structure water.pdb
    number {waternum}
    inside box {waterloc[0]} {waterloc[1]} {waterloc[2]} {waterloc[3]} {waterloc[4]} {waterloc[5]}
end structure

{stru2}
    """

    output1 = ''
    for i in range(no_surfts):
        struout = template1.format(
            surfs_name = surfts_type[i],
            num_surfs = itf_surfts[i],
            surfloc = surflocu[i],
            hydrophilic = hindexlst[i],
            hydrophilic_anchor = anchorlocu[i][0],
            hydrophobic = tindexlst[i],
            hydrophobic_anchor = anchorlocu[i][1],
        )
        output1 += struout


    output2 = ''
    for i in range(no_surfts):
        struout = template2.format(
            surfs_name = surfts_type[i],
            num_surfs = itf_surfts[i],
            surfloc = surflocd[i],
            hydrophilic = hindexlst[i],
            hydrophilic_anchor = anchorlocd[i][0],
            hydrophobic = tindexlst[i],
            hydrophobic_anchor = anchorlocd[i][1],
        )
        output2 += struout


    packmolinp = templateT.format(
        output_filename = simulationsys_name,
        waternum = now,
        waterloc = waterloc,
        stru1 = output1,
        stru2 = output2
    )
    packinpfn = simulationsys_name+'.inp'
    with open(packinpfn, 'w') as output_file:
        output_file.write(packmolinp)

def ff_unilier(surfts_type,simulationsys_name):
    templateff = """;Genarate by MFST 
[ atomtypes ]
{atomtypes}

{otherinfo}
"""

    atminp = ''
    othinp = ''

    for surfts in surfts_type:
        ff = surfts_ff_dict[surfts]
        with open (ff,'r') as file:
            lines = file.read()
        newlines = lines.replace('opls',surfts)
        atomtypes_match = re.search(r'\[ atomtypes \](.*?)\[', newlines, re.DOTALL)
        atomtypes_content = atomtypes_match.group(1).strip() if atomtypes_match else None
        # 使用正则表达式选择 [moleculetype] 及其以下的全部内容
        moleculetype_match = re.search(r'(\[ moleculetype \].*?)$', newlines, re.DOTALL)

        # 获取匹配到的内容
        moleculetype_content = moleculetype_match.group(1).strip() if moleculetype_match else None
        atminp += atomtypes_content
        othinp += moleculetype_content

    ffunit = templateff.format(
        atomtypes = atminp,
        otherinfo = othinp
    )
    fffn = simulationsys_name+'.itp'
    with open(fffn, 'w') as output_file:
        output_file.write(ffunit)

def move_file_to_directory(file, new_directory):
    # 创建目录或判断目录是否存在
    if os.path.exists(new_directory):
        user_input = input(f"The directory '{new_directory}' already exists. Do you want to move the file? (y/n): ").lower()
        if user_input == 'y':
            n = 1
        elif user_input == 'n':
            # 建立带数字编号的新目录
            index = 1
            while os.path.exists(new_directory):
                new_directory = f'{new_directory}_{index}'
                index += 1
            os.makedirs(new_directory)
    else:
        os.makedirs(new_directory)
    # 移动文件到新目录
    shutil.move(file, new_directory)

    print(f"File moved to {new_directory}")

def run_packmol(packmolpath,input_file,filepath):
    packmol_path = packmolpath  # 替换为你的Packmol可执行文件路径
    directory, file_name = os.path.split(input_file)
    os.chdir(directory)
    # 构建Packmol命令
    packmol_cmd = [packmol_path, '<', file_name]

    with open(file_name, 'r') as f:
        subprocess.run([packmol_path], stdin=f, shell=True)
        
    os.chdir(filepath)

def generate_topfile(simulationsys_name,WFF,surfts_type,itf_surfts,now):
    templateTOP="""#include "oplsaa.ff/forcefield.itp"
#include "{simulationsys}.itp"
#include "oplsaa.ff/{WFF}.itp"
#include "oplsaa.ff/ions.itp"


[ system ]
;name
{simulationsys} Simulationbox

[ molecules ]
{SurftsInfo}
SOL      {WatNum}
{SurftsInfo}
"""
    SurftsInfo = ''
    for surft in range(len(surfts_type)):
        singeln = ''
        singeln += surfts_type[surft]
        singeln += '\t'
        singeln += str(itf_surfts[surft])
        singeln += '\n'
        SurftsInfo += singeln
    SurftsInfo = SurftsInfo.rstrip('\n')
    outputTOP = templateTOP.format(
        simulationsys = simulationsys_name,
        WFF = WFF,
        SurftsInfo = SurftsInfo,
        WatNum = now,
    )
    topfn = simulationsys_name+'.top'
    with open(topfn, 'w') as output_file:
        output_file.write(outputTOP)

def pyModeling_main():
    logging.basicConfig(filename='pyModeling.log', level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    start_time = time.time()

    config_file_path = './MFST.rc'  
    no_surfts, surfts_type, interfc_type,itf_surfts, waterbox_size_list, simulationbox_list, simulationsys_name,temp,packmolpath,WFF = read_config(config_file_path)


    # 检查 surfts_type 的变量个数与 no_surfts 是否一致
    if len(surfts_type) != no_surfts:
        print("Error: Number of surfts_type variables does not match no_surfts. Program will terminate.")
        exit(1)

    # 配置映射
    load_surfts_type_mapping()

    # 分析亲水头基和疏水尾基
    sfpath = []
    ffpath = []
    for surfts in surfts_type:
        sfpath.append(surfts_type_dict[surfts])
        ffpath.append(surfts_ff_dict[surfts])


    hindexlst = []
    tindexlst = []
    distancelst = []
    for sf in sfpath:
        hindex, tindex, distance = find_hydrophobic_tail_atoms(sf)
        hindexlst.append(hindex)
        tindexlst.append(tindex)
        distancelst.append(distance)

    if os.path.exists(simulationsys_name):
        n =1
    else:
        os.makedirs(simulationsys_name)

    for i in range(len(sfpath)):
        newname = surfts_type[i]+'.pdb'
        destination_path = os.path.join(simulationsys_name, newname)
        shutil.copy2(sfpath[i], destination_path)
    shutil.copy2('./database/sysmol/water.pdb', simulationsys_name)

    # 计算水分子个数
    now = get_number_of_water(waterbox_size_list,temp)


    generate_packmol_input(no_surfts, surfts_type, interfc_type,itf_surfts, waterbox_size_list, simulationbox_list, simulationsys_name,hindexlst,tindexlst,distancelst,now)
    ff_unilier(surfts_type,simulationsys_name)
    generate_topfile(simulationsys_name,WFF,surfts_type,itf_surfts,now)

    move_file_to_directory(simulationsys_name+'.inp',simulationsys_name)
    move_file_to_directory(simulationsys_name+'.itp',simulationsys_name)
    move_file_to_directory(simulationsys_name+'.top',simulationsys_name)
    script_path = os.path.abspath(sys.argv[0])
    script_directory = os.path.dirname(script_path)
    run_packmol('/home/Data/zhaoxz/application/packmol/packmol/packmol',simulationsys_name +'/'+simulationsys_name+'.inp',script_directory)

    
    absolute_path = os.path.abspath(simulationsys_name)
    workdir = absolute_path
    logging.info(f'Project Name: {simulationsys_name}')
    logging.info(f'Working Directory: {workdir}')
    end_time = time.time()
    elapsed_time = end_time - start_time
    logging.info(f'Program finished in {elapsed_time:.2f} seconds.')
    logging.info('MFST wishes you: Have a great day! :-) ')
        
if __name__ == "__main__":
    pyModeling_main()