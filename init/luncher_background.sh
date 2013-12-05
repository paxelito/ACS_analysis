echo "--> Simulations archive extraction..."
tar -xvzf exps.tar.gz

echo "--> CaRNeSS compiling..."
cd ./exps/sourceCode/
g++ -Wall -ansi -lm -o ../SIMS/carness *.cpp
cd ../SIMS/
chmod 755 carness

echo "--> Simulations..."
chmod 755 *simulation.sh
./1_simulation.sh &
./2_simulation.sh &
#./3_simulation.sh &
#./4_simulation.sh &
#./5_simulation.sh &
#./6_simulation.sh &
#./7_simulation.sh &
#./8_simulation.sh &
