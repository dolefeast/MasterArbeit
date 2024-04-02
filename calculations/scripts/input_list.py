def input_list(data_type=int):
    """
    Asks user for inputs for a list.
    Parameters:
        data_type: the data type to force the data to
    """
    array = []
    item = input("Insert data (empty input to stop)...")
    while item != "":
        try:
            array.append(data_type(item))
        except ValueError:
            print(f"{item} is not {data_type}!, try again.")
        item = input("Insert data...")

    return array
