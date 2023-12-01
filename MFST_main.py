# main_program.py

import pyModeling
import pySimulating
import pyAnalysing

def MSmain():
    # 调用 pyModeling 模块
    print("Running pyModeling module...")
    pyModeling.pyModeling_main()

    # 调用 pySimulating 模块
    print("Running pySimulating module...")
    pySimulating.pySimulating_main()

def Amain():
    # 调用 pyAnalysing 模块
    print("Running pyAnalysing module...")
    pyAnalysing.pyAnalysing_main()

if __name__ == "__main__":
    print("""
 _  _  ____  ____  ____ 
( \/ )(  __)/ ___)(_  _)
/ \/ \ ) _) \___ \  )(  
\_)(_/(__)  (____/ (__) 
""")
    print('MFST MolFilmStabTool © 2023 MTSD@UPC')
    option = input('Which part you want? MS or A?')
    if option == 'MS':
        MSmain()
    elif option == 'A':
        Amain()
    else:
        print('Sorry! You can only type MS or A')
