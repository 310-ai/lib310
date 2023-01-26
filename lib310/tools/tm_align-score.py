import os
import subprocess



def tm_allign( pdb_file_1,pdb_file_2):
    import os.path
    if (os.path.exists(r'TMalign')) == False:
        os.system('g++ -O3 -ffast-math -lm -o TMalign TMalign.cpp')
    batcmd = './TMalign {file1} {file2}'.format(file1=pdb_file_1, file2=pdb_file_2)
    # result = subprocess.check_output(batcmd, shell=True)
    result2 = subprocess.getoutput(batcmd)
    # print("result::: ", result)
    # result2[result2.find('TM-score    '):result2.find('TM-score    ') + 32]
    return result2



def tm_score( pdb_file_1,pdb_file_2):
    import os.path
    if (os.path.exists(r'TMscore')) == False:
        os.system('g++ -O3 -ffast-math -lm -o TMscore TMscore.cpp')
    batcmd = './TMscore {file1} {file2}'.format(file1=pdb_file_1, file2=pdb_file_2)
    # result = subprocess.check_output(batcmd, shell=True)
    result2 = subprocess.getoutput(batcmd)
    # print("result::: ", result)
    # result2[result2.find('TM-score    '):result2.find('TM-score    ') + 32]
    return result2




