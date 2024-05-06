def float_to_str(value, sig_digs):
    return str(
        round(
            float(
                value
            ),
            sig_digs
        )
    ).replace(".", "_")

def str_to_float(value:str):
    return float(
        value.replace("_", ".")
        )
