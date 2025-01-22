from pathlib import Path, PurePath
import re

from minimal_working_example import str_to_float


bcs = "dirichlet"
quantity_names = ['eigenvalues', 'eigenstates', 'eigenstate_gradients', 'rho', 'A0_induced']
def get_directory_m_a(bcs, read_mode='str'):

    dirs = Path("saved_solutions")

    [print(f"[{i+1}]", dir_name.name) for i, dir_name in enumerate(dirs.glob("*"))]

    while True:
        try:
            index = int(input("Choose the desired directory... "))
            directory = list(dirs.glob('*'))[index-1]
            break
        except ValueError:
            print("Input must be a natural number, try again")
        except IndexError:
            print(f"Input must be a number from 1 to {len(list(dirs.glob('*')))}")


    for name in ["rho"]: # Just need one of those
        quantity_dir = Path(PurePath(directory) / bcs / name)

        m_a_set = set()
        for filename in quantity_dir.glob("**/*txt"):
            try:
                if read_mode=='str':
                    matches = re.findall("\d+_\d+", str(filename.name))
                    m, a, lambda_value = matches
                elif read_mode=='float':
                    m, a, lambda_value = map(str_to_float, re.findall("\d+_\d+", str(filename.name)))
                else:
                    raise ValueError("read_mode must be either 'str', 'float' but it was " + read_mode)
            except ValueError as e:
                print(e)
                continue
            m_a_set.add((m,a))

    if len(m_a_set) > 1:
        print("The distinct m, a values are:")
        for i, m_a in enumerate(m_a_set):
            print(f"[{i+1}]: m = {list(m_a)[0]}, a = {list(m_a)[1]}") 
        while True:
            try:
                index = int(input("Choose the desired m, a value ... "))
                m_a_choice = list(m_a_set)[index-1]
                break
            except ValueError:
                print("Input must be a natural number, try again")
            except IndexError:
                print(f"Input must be a number from 1 to {len(list(dirs.glob('*')))}")
    else:
        print(f"There is only one distinct case (m, a) = {list(m_a_set)[0]}")
        m_a_choice = list(m_a_set)[0]

    return directory.name, m_a_choice[0], m_a_choice[1]

if __name__ == "__main__":
    bcs = "dirichlet"
    print(get_directory_name(bcs))
