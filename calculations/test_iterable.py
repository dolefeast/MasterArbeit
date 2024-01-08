from math_objects.iterable_output import iterable_output

@iterable_output
def identity(x):
    return x

print(list(map(identity, range(5))))
