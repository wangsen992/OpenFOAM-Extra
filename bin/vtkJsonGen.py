#!/bin/python3

import json
import os
import pathlib

def organize_vtk(postProcPath):
    print("Collecting postProcessing files into json format")

    if type(postProcPath) != pathlib.Path:
        postProcPath = pathlib.Path(postProcPath)
    resultDict = {}
    for fpath in postProcPath.rglob("*.vtk"):
        print(fpath)
        try:
            resultDict[fpath.stem].append(dict(name=fpath.as_posix(),
                time=float(fpath.parent.name)))
        except KeyError:
            resultDict[fpath.stem] = list()
            resultDict[fpath.stem].append(dict(name=fpath.as_posix(),
                time=float(fpath.parent.name)))
    
    return resultDict

if __name__ == "__main__":

    resultDict = organize_vtk("./postProcessing")
    
    for key in resultDict.keys():
        print(f"Writing series file for {key}")
        with open(key + ".vtk.series", 'w') as f:
            keyDict = {"file-series-version": "1.0", 
                              "files" : resultDict[key]}
            json.dump(keyDict, f)




