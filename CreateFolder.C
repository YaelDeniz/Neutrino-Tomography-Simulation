#include <iostream>
#include <stdio.h>
#include <sys/stat.h>

void CreateFolder(){
std::string ResultsPATH = "/home/dehy0499/NuOscillation-Tomography/Neutrino-Tomography-Simulation/";
std::string DirName = ResultsPATH+"Test777";
std::string Folder(DirName);
const char* ccx = Folder.c_str();

mkdir(ccx, 0777);
}
