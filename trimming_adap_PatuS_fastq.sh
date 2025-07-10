#!/bin/bash

# Ruta de entrada: carpeta con los archivos de secuencias
input_dir="../01.QC/PatuS_reads"

# Ruta de salida: carpeta donde guardarás los archivos procesados
output_dir="../01.QC/PatuS_reads/adapt_trimming"
mkdir -p "$output_dir"  # Crear la carpeta de salida si no existe

# Archivo de adaptadores en formato FASTA
adapter_file="../TFM_previo/adapters.fa"

# Procesar todos los archivos pareados (_R1_ y _R2_)
for r1 in "$input_dir"/*_R1_*.fastq.gz; do
    # Encontrar el archivo _R2_ correspondiente
    r2=${r1/_R1_/_R2_}

    # Verificar que los archivos existan
    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "Error: Uno de los archivos (_R1_ o _R2_) no existe. Saltando..."
        echo "Archivo forward (R1): $r1"
        echo "Archivo reverse (R2): $r2"
        continue
    fi

    # Obtener el nombre base del archivo
    base_name_r1=$(basename "$r1")
    base_name_r2=$(basename "$r2")

    echo "Procesando $base_name_r1 y $base_name_r2 con fastp..."

    # Ejecutar fastp para recortar adaptadores, usando los mismos nombres de entrada para la salida
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$output_dir/$base_name_r1" \
        -O "$output_dir/$base_name_r2" \
        --adapter_fasta "$adapter_file" \
        --disable_quality_filtering \
        --disable_length_filtering

    echo "$base_name_r1 y $base_name_r2 han sido procesados."
done

echo "Todos los archivos han sido procesados y guardados en $output_dir."
