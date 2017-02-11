cd src
touch ./*
make
cd ..
rm output/*
./bin/OneDimWaveEqn
./Animate output/*