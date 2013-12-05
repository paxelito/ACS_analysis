echo "--> Simulations archive extraction..."
tar -xvzf exps.tar.gz

echo "--> CaRNeSS compiling..."
cd ./exps/sourceCode/
g++ -Wall -ansi -lm -o ../SIMS/carness *.cpp
cd ../SIMS/
chmod 755 carness

echo "--> Simulations..."
chmod 755 *simulation.sh
screen -d -m -S C1 ./1_simulation.sh
screen -d -m -S C2 ./2_simulation.sh 
screen -d -m -S C3 ./3_simulation.sh 
screen -d -m -S C4 ./4_simulation.sh
screen -d -m -S C5 ./5_simulation.sh 
screen -d -m -S C6 ./6_simulation.sh 
screen -d -m -S C7 ./7_simulation.sh 
screen -d -m -S C8 ./8_simulation.sh 
