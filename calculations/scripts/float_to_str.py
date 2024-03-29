def float_to_str(value, sig_digs):
    return str(round(float(value), sig_digs)).replace(".", "_")
