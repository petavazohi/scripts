def get_type(item):
    if isinstance(item, dict):
        ret = {}
        for key in item:
            ret[key] = get_type(item[key])
    elif isinstance(item, list) or isinstance(item, tuple):
        ret = []
        for x in item:
            ret.append(get_type(x))
        if isinstance(item, tuple):
            ret = tuple(ret)
    else:
        ret = type(item)
    return ret
