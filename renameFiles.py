import os, re
from glob import glob
direc = r'Z:\visitor\ch6617\bm31\20230329/'

def run(direc):
    os.chdir(direc)
    digits = 3
    for root, dir, files in os.walk(os.getcwd()):
        if 'columns' in root and not 'merge' in root and 'norm' not in root:

            os.chdir(root)
            print(root)
            datfiles = glob('*.dat')
            norfiles = glob('*.nor')
            allfiles = datfiles + norfiles
            for file in allfiles:

                fileno = re.split('\.|_',file)[-2]
                intfileno = int(fileno)
                newfilename = file.replace(f'_{fileno}.',f'_{intfileno:0{digits}d}.')

                try:
                    os.rename(file,newfilename)
                except FileExistsError:
                    os.remove(file)
                    print(f'{root}/{file}')
        elif 'columns' in root and 'norm' in root:
            os.chdir(root)
            ffiles = glob('*F.nor')
            tfiles = glob('*T.nor')
            for file in ffiles:
                fileno = re.split('\.|_',file)[-2].replace('F','')
                intfileno = int(fileno)
                newfilename = file.replace(f'_{fileno}F.',f'_{intfileno:0{digits}d}F.')
                try:
                    os.rename(file,newfilename)
                except FileExistsError:
                    os.remove(file)
                    print(f'{root}/{file}')
            for file in tfiles:
                fileno = re.split('\.|_',file)[-2].replace('T','')
                intfileno = int(fileno)
                newfilename = file.replace(f'_{fileno}T.',f'_{intfileno:0{digits}d}T.')
                try:
                    os.rename(file,newfilename)
                except FileExistsError:
                    os.remove(file)
                    print(f'{root}/{file}')
if __name__ == '__main__':
    run(direc = direc)
