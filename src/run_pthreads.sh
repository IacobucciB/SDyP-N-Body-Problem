#!/bin/bash

echo "Seleccione el sistema operativo para compilar:"
select os_choice in "Linux" "Windows"; do
    case $os_choice in
        Linux)
            gcc -o pthreads pthreads.c -pthread -lm
            break
            ;;
        Windows)
            gcc -o pthreads pthreads.c -lpthread -lm
            break
            ;;
        *)
            echo "Opción inválida. Intente de nuevo."
            ;;
    esac
done

if [ $? -ne 0 ]; then
    echo "Error: Compilation failed."
    exit 1
fi

./pthreads 512 200 1000 4 > output_pthreads.txt
if [ $? -ne 0 ]; then
    echo "Error: Execution failed."
    exit 1
fi

echo "Execution completed successfully."