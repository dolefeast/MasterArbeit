def get_lambda_index(lambda_array, lambda_target):
    if lambda_target > max(lambda_array):
        return max(lambda_array)
    if lambda_target in lambda_array:
        return lambda_array.index(lambda_target)
    else:
        for i, lambda_value in enumerate(lambda_array):
            if lambda_value>lambda_target:
                return i-1
