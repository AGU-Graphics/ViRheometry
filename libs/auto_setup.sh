brew install boost-python3
brew install eigen
python3 -m pip install numpy
python3 -m pip install taichi
brew install cmake
cd libs
git clone https://github.com/CMA-ES/libcmaes
cp lcmaes.cc libcmaes/python/lcmaes.cc
cd libcmaes
mkdir build
cd build
cmake -DLIBCMAES_BUILD_PYTHON=ON ..
make -j8
cp ./python/lcmaes.so ../../../libs
cd ..
cd ..
cd ..
python3 -m pip install numpy 
python3 -m pip install matplotlib
python3 -m pip install plotly
python3 -m pip install tqdm
python3 -m pip install numdifftools
python3 -m pip install pandas
python3 -m pip install streamlit
brew install tcl-tk
brew install opencv
brew install python-tk
