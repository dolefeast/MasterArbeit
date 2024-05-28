def main(
        m,
        lambda_value,
        ):
    from physics import Vacuum_Polarization

    import matplotlib.pyplot as plt
    from multiprocessing import Process

    fig, ax = plt.subplots(figsize=(16, 9))

    system = Vacuum_Polarization(
            lambda_value=lambda_value,
            m=m,
            )

    for i in range(3):
        system.update_eigenstates()
        ax.plot(
                system.z,
                system.rho,
                alpha = (1+i)/3,
                )
    plt.show()

if __name__ == "__main__":
    from multiprocessing import Process

    def bubble_sort(array):
        for _ in range(int(5e6)):
            check = True
            while check == True:
              check = False
              for i in range(0, len(array)-1):
                if array[i] > array[i+1]:
                  check = True
                  temp = array[i]
                  array[i] = array[i+1]
                  array[i+1] = temp
            print("Array sorted: ", array)

    p = Process(target=bubble_sort, args=([1,9,4,5,2,6,8,4],))
    p.start()
    p.join()
