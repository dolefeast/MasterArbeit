import regex as re 

from utils.str_to_float import str_to_float 
float_re = re.compile("\d+_\d+")

def get_lambda_value(filename):
    # Each filename has only the mass and lambda_value parameter
    return str_to_float(float_re.findall(str(filename))[1])
