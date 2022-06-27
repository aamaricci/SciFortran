#Building scifor
#Errors
set -e

cd scifor
mkdir build
cd build

echo "cmake .."
cmake ..

echo "make"
make

echo "make install"
make install

echo "source ~/opt/scifor/gnu/*/bin/scifor_config_user.sh" >> ~/.scifor_config_user
echo -e "\e[32m scifor installed and sourced \e[0m"

