for i in 10; 
do
    echo ${i}
    mkdir -p K${i}
    cd K${i}
    mkdir -p vecs
    screen -d -m ../porous ${i} 1 10
    cd ..
done
