for i in 0.1 0.2 0.3 0.5; 
do
    echo ${i}
    mkdir -p K=${i}
    cd K=${i}
    mkdir -p vecs
    screen -d -m ../porous ${i}
    cd ..
done
