cd src
touch ./*
make
cd ..
rm output/*
./bin/TwoDimWaveEqn
./Animate output/*
