To install app, navigate to for_redistribution file and run MyAppInstaller_web.install
If you wish to install in /usr/local (default) must run as sudo user
sudo ./MyAppInstaller_web.install

To add matlab to the path add this line of code at end of ~/.bashrc script
 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/MATLAB/MATLAB_Runtime/v90/bin/glnxa64

where /usr/local/MATLAB/ is the destination of the installation

direct questions to christian.t.meyer@vanderbilt.edu


