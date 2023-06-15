import json

def read_params(params_json: str) -> dict:
    params = json.load(open(params_json))
    return params