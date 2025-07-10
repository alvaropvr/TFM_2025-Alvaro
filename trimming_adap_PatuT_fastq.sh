#!/bin/bash

# Ruta de entrada: carpeta con los archivos de secuencias
input_dir="../01.QC/PatuT_reads"

# Ruta de salida: carpeta donde guardarás los archivos procesados
output_dir="../01.QC/PatuT_reads/trimming_adap"
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

    # Obtener el nombre base del archivo sin extensiones
    base_name=$(basename "$r1" _R1_001.fastq.gz)

    echo "Procesando $base_name con fastp..."

    # Ejecutar fastp para recortar adaptadores
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$output_dir/${base_name}_R1_trimmed.fastq.gz" \
        -O "$output_dir/${base_name}_R2_trimmed.fastq.gz" \
        --adapter_fasta "$adapter_file" \
        --disable_quality_filtering \
        --disable_length_filtering

    echo "$base_name ha sido procesado."
done

echo "Todos los archivos han sido procesados y guardados en $output_dir."
