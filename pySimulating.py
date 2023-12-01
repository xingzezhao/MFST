#!/usr/bin/env python
# coding: utf-8

import configparser
import os
import subprocess

def convert_to_absolute_path(path,MFST_path):
    if not os.path.isabs(path):
        current_dir = MFST_path
        path = os.path.join(current_dir, path)
    return path

def read_config(config_file):
    # 创建配置解析器对象
    config = configparser.ConfigParser()
    MFST_path = os.path.abspath('.')
#     print(MFST_path)
    # 读取配置文件
    config.read(config_file)
    for section in config.sections():
        for option, value in config.items(section):
            if value[0] == '.':
                config.set(section, option, os.path.abspath(value))
    # 获取配置项的值
    config_values = {
        'no_surfts': config.getint('pyModeling input', 'no_surfts'),
        'surfts_type': config.get('pyModeling input', 'surfts_type').split(),
        'interfc_type': config.get('pyModeling input', 'interfc_type'),
        'itf_surfts': [int(value) for value in config.get('pyModeling input', 'itf_surfts').split()],
        'waterbox_size_list': [float(value) for value in config.get('pyModeling input', 'waterbox_size').split()],
        'simulationbox_list': [float(value) for value in config.get('pyModeling input', 'simulationbox_size').split()],
        'simulationsys_name': config.get('pyModeling input', 'simulationsys_name'),
        'temp': config.getfloat('pySimulating input', 'Temp'),
        'packmolpath': config.get('pyModeling input', 'packmolpath'),
        'WFF': config.get('pySimulating input', 'SimuWatMod'),
        'eqnsteps':  config.getint('pySimulating input', 'eqnsteps'),
        'simnsteps':  config.getint('pySimulating input', 'simnsteps'),
        'simdt': config.getfloat('pySimulating input', 'dt'),
        'outflames': config.getint('pySimulating input', 'outflames'),
        'rcoulomb': config.getfloat('pySimulating input', 'rcoulomb'),
        'rlist': config.getfloat('pySimulating input', 'rlist'),
        'rvdw': config.getfloat('pySimulating input', 'rvdw'),
        'gmxpath': config.get('pySimulating input', 'gmxpath'),
        'pions': config.get('pySimulating input', 'pions'),
        'nions': config.get('pySimulating input', 'nions'),
        'conc': config.getfloat('pySimulating input', 'conc'),
        'emtol': config.get('pySimulating input', 'emtol'),
        'jcrls': config.get('pySimulating input', 'jcrls'),
    }
    
    jcrtls = config_values['jcrls']
    if jcrtls == 'Slurm':
        config_values['jcspath'] = config.get('pySimulating input', 'sltmpath')
        config_values['jcline'] = 'sbatch'
    elif jcrtls == 'PBS':
        config_values['jcspath'] = config.get('pySimulating input', 'pbstmpath')
        config_values['jcline'] = 'qsub'
    else:
        print(f'ERROR: Cannot recognize the Job Crtl Sys:{jcrtls}, you should use Slurm or PBS')

    return config_values

def generate_mdpfile(config_values):
    ## 编辑离子化和能量最小化参数文件
    emtol = config_values['emtol']
    rcoulomb = config_values['rcoulomb']
    rvdw = config_values['rvdw']
    eqnsteps = config_values['eqnsteps']
    simnsteps = config_values['simnsteps']
    dt = config_values['simdt']
    nst = config_values['nst']
    nstlog = config_values['nstlog']
    rlist = config_values['rlist']
    rcoulomb = config_values['rcoulomb']
    rvdw = config_values['rvdw']
    temp = config_values['temp']
    simulationsys_name = config_values['simulationsys_name']


    templateEM = """integrator  	= steep
emtol       	= {emtol}
emstep      	= 0.01
nsteps      	= 5000
nstlist         = 1
cutoff-scheme	= Verlet
ns_type         = grid
coulombtype     = cutoff
rcoulomb        = {rcoulomb}
rvdw            = {rvdw}
pbc             = xyz
    """
    outputEM = templateEM.format(
        emtol = emtol,
        rcoulomb = rcoulomb,
        rvdw = rvdw,
    )
    emfn = simulationsys_name+'_em.mdp'
    with open(emfn, 'w') as output_file:
        output_file.write(outputEM)
        
    templateEQ = """title 					= NVT equilibration for system
integrator				= md
nsteps					= {eqnsteps}
dt		    			= {dt}
comm-mode				= Linear
comm-grps				= system
nstxout					= {nst}
nstvout					= {nst}
nstenergy				= {nst}
nstlog					= {nstlog}
energygrps				= system
xtc-grps				= system
continuation			= no
constraint_algorithm	= lincs
constraints				= H-bonds
lincs_iter				= 1
lincs_order				= 4
ns_type					= grid
nstlist					= 200
rlist					= {rlist}
rcoulomb				= {rcoulomb}
rvdw					= {rvdw}
cutoff-scheme			= Verlet
coulombtype				= PME
pme_order				= 4
fourierspacing			= 0.16
tcoupl					= V-rescale
tc-grps					= system
tau_t					= 0.1
ref_t					= {temp}
pcoupl					= no
pbc						= xyz
DispCorr				= EnerPres
gen_vel					= yes
gen_temp				= {temp}
gen_seed				= 9959
    """
    outputEQ = templateEQ.format(
        eqnsteps = eqnsteps,
        dt = dt,
        nst = nst,
        nstlog = nstlog,
        rlist = rlist,
        rcoulomb = rcoulomb,
        rvdw = rvdw,
        temp = temp,
    )
    eqfn = simulationsys_name+'_eq.mdp'
    with open(eqfn, 'w') as output_file:
        output_file.write(outputEQ)
        
    templateMD = """title 					= NVT equilibration for system
integrator				= md
nsteps					= {simnsteps}
dt		    			= {dt}
comm-mode				= Linear
comm-grps				= system
nstxout					= {nst}
nstvout					= {nst}
nstenergy				= {nst}
nstlog					= {nstlog}
energygrps				= system
xtc-grps				= system
continuation			= yes
constraint_algorithm	= lincs
constraints				= H-bonds
lincs_iter				= 1
lincs_order				= 4
ns_type					= grid
nstlist					= 200
rlist					= {rlist}
rcoulomb				= {rcoulomb}
rvdw					= {rvdw}
cutoff-scheme			= Verlet
coulombtype				= PME
pme_order				= 4
fourierspacing			= 0.16
tcoupl					= V-rescale
tc-grps					= system
tau_t					= 0.1
ref_t					= {temp}
pcoupl					= no
pbc						= xyz
DispCorr				= EnerPres
gen_vel					= no
    """
    outputMD = templateMD.format(
        simnsteps = simnsteps,
        dt = dt,
        nst = nst,
        nstlog = nstlog,
        rlist = rlist,
        rcoulomb = rcoulomb,
        rvdw = rvdw,
        temp = temp,
    )
    mdfn = simulationsys_name+'_md.mdp'
    with open(mdfn, 'w') as output_file:
        output_file.write(outputMD)

def find_workdir(pymdlog,simulationsys_name):
    with open(pymdlog, 'r') as log_file:
        # 遍历文件的每一行
        for line in log_file:
            if 'Project Name' in line:
                if line.split(': ')[-1].strip() == simulationsys_name :
                    # 检查是否包含目录信息
                    working_directory_line = next(log_file).strip()
                    if 'Working Directory' in working_directory_line:
                        # 提取目录路径
                        directory_path = working_directory_line.split(': ')[-1].strip()

                        # 切换工作目录
                        os.chdir(directory_path)
                        print(f'Changed working directory to: {directory_path}')
                        break  # 读到目录信息后退出循环
                
def pdb2box_ion_em(config_values):
    gmx_path = config_values['gmxpath']
    inistr = config_values['simulationsys_name'] + '.pdb'
    inibox = config_values['simulationsys_name'] + 'box.gro'
    systop = config_values['simulationsys_name'] + '.top'
    gmx_cmd1 = [gmx_path,'editconf' ,'-f', inistr, '-o', inibox, '-box', str(config_values['simulationbox_list'][0]), str(config_values['simulationbox_list'][1]), str(config_values['simulationbox_list'][2])]
    emmdppath = config_values['simulationsys_name'] + '_em.mdp'
    eqmdppath = config_values['simulationsys_name'] + '_eq.mdp'
    mdmdppath = config_values['simulationsys_name'] + '_md.mdp'
    iontpr = config_values['simulationsys_name'] + '_ions.tpr'
    gmx_cmd2 = [gmx_path, 'grompp', '-f', emmdppath, '-c', inibox,
               '-o', iontpr,'-p',systop, '-maxwarn','1']
    iongro = config_values['simulationsys_name'] + '_ions.gro'
    gmx_cmd3 = [gmx_path, 'genion', '-s', iontpr, '-p', systop, '-o', iongro,
              '-pname', config_values['pions'],  '-nname', config_values['nions'],]
    if config_values['conc'] == 0:
        gmx_cmd3.append('-neutral')
    else:
        gmx_cmd3.append('-conc')
        gmx_cmd3.append(config_values['conc'])

    emtpr = config_values['simulationsys_name'] + '_em.tpr'
    eminp = config_values['simulationsys_name'] + '_em'
    gmx_cmd4 = [gmx_path, 'grompp', '-f', emmdppath, '-c', iongro,
               '-o', emtpr,'-p',systop]
    gmx_cmd5 = [gmx_path, 'mdrun', '-deffnm', eminp]
    subprocess.run(gmx_cmd1)
    subprocess.run(gmx_cmd2)
    process = subprocess.Popen(gmx_cmd3, stdin=subprocess.PIPE, text=True)
    process.communicate(input='SOL\n')
    process.wait()
    subprocess.run(gmx_cmd4) 
    subprocess.run(gmx_cmd5)    

def run_nvtnmd(config_values):
    
    jcline = config_values['jcline']
    jcspath = config_values['jcspath']
    with open(jcspath) as file:
        jcrlfl = file.read()

    gmx_path = config_values['gmxpath']
    ordline = f'source {gmx_path}'
    ordline+= '\n'
    systop = config_values['simulationsys_name'] + '.top'
    emgro = config_values['simulationsys_name'] + '_em.gro'
    eqmdp = config_values['simulationsys_name'] + '_eq.mdp'
    eqtpr = config_values['simulationsys_name'] + '_eq.tpr'
    eqol = config_values['simulationsys_name'] + '_eq'
    ordline+= f'gmx grompp -f {eqmdp} -c {emgro} -r {emgro} -o {eqtpr} -p {systop}'
    ordline+= '\n'
    ordline+= f'gmx mdrun -deffnm {eqol}'
    eqgro = config_values['simulationsys_name'] + '_eq.gro'
    mdmdp = config_values['simulationsys_name'] + '_md.mdp'
    mdtpr = config_values['simulationsys_name'] + '_md.tpr'
    mdol = config_values['simulationsys_name'] + '_md'
    ordline+= '\n'
    ordline+= f'gmx grompp -f {mdmdp} -c {eqgro} -r {eqgro} -o {mdtpr} -p {systop}'
    ordline+= '\n'
    ordline+= f'gmx mdrun -deffnm {mdol}'

    jcrlflop = jcrlfl.format(
        orderline = ordline
    )

    mdfn = config_values['simulationsys_name']+'.slurm'
    with open(mdfn, 'w') as output_file:
        output_file.write(jcrlflop)

    jcsubln = [jcline,mdfn]
    subprocess.run(jcsubln) 

def pySimulating_main():
    config_file_path = './MFST.rc'  
    config_values = read_config(config_file_path)
    MFST_path = os.path.abspath('.')
    find_workdir('./pyModeling.log',config_values['simulationsys_name'])


    nst = round(config_values['simnsteps'] / config_values['outflames'])
    nstlog = int(2*nst)
    config_values['nst'] = nst
    config_values['nstlog'] = nstlog
    generate_mdpfile(config_values)
    pdb2box_ion_em(config_values)
    run_nvtnmd(config_values)

    os.chdir(MFST_path)

if __name__ == "__main__":
    pySimulating_main()
