def floatToStr(self, value):
    """
    Changes the float 1.234 to "1_234"
    """
    if isinstance(value, str):
        return value
    return str(
        round( float( value), self.sigDigs)
    ).replace(".", "_")

def strToFloat(self, value:str):
    return float(
        value.replace("_", ".")
        )
