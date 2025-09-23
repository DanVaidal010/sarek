#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import math
import argparse
import traceback

def obtener_insert_metrics(json_path):
    """
    Calcula la media y desviación estándar del tamaño de inserto
    a partir del histograma presente en un archivo `.fastp.json`.
    """
    try:
        # Cargar archivo JSON de fastp
        with open(json_path) as f:
            data = json.load(f)

        # Extraer histograma de tamaños de inserto
        insert_data = data.get("insert_size", {})
        hist = insert_data.get("histogram", [])

        # Desanidar si es una lista de listas
        if hist and isinstance(hist[0], list):
            hist = [item for sublist in hist for item in sublist]

        # Calcular número total de lecturas
        total_reads = sum(hist)
        if total_reads == 0:
            raise ValueError("Histograma vacío o no válido.")

        # Media ponderada: sum(i * count) / total
        media = sum(i * count for i, count in enumerate(hist)) / total_reads

        # Varianza ponderada y desviación estándar
        variance = sum(((i - media) ** 2) * count for i, count in enumerate(hist)) / total_reads
        std_dev = round(math.sqrt(variance), 2)

        return {
            "Mean_insert": round(media, 2),
            "Insert_SD": std_dev
        }

    except Exception as e:
        # En caso de error, registrar y devolver valores nulos
        print(f"❌ Error al procesar archivo JSON: {e}")
        traceback.print_exc()
        return {
            "Mean_insert": None,
            "Insert_SD": None
        }

if __name__ == "__main__":
    # Argumentos CLI
    parser = argparse.ArgumentParser(description="Calcula media y SD del tamaño de inserto desde .fastp.json")
    parser.add_argument("-j", "--json", required=True, help="Ruta al archivo .fastp.json")
    parser.add_argument("-o", "--output", required=True, help="Archivo de salida (ej: sample.insert.txt)")
    args = parser.parse_args()

    # Calcular métricas
    resultados = obtener_insert_metrics(args.json)

    # Guardar en archivo de texto
    with open(args.output, 'w') as out:
        for clave, valor in resultados.items():
            out.write(f"{clave}\t{valor}\n")

    print(f"✅ Inserto (media y SD) guardado en {args.output}")
