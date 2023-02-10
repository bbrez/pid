#!/usr/bin/env python3

from PIL import Image
import numpy
import argparse
import heapq
import struct

img_width = 0
img_height = 0
is_gray = False

# Função que calcula a frequência das cores
def calc_freq(array):
    freq = {}
    for i in array:
        if i in freq:
            freq[i] += 1
        else:
            freq[i] = 1
    return freq

# Função que constroi a árvore de Huffman
def build_tree(freqs):
    heap = [[wt, [sym, '']] for sym, wt in freqs.items()]
    heapq.heapify(heap)
    while len(heap) > 1:
        lo = heapq.heappop(heap)
        hi = heapq.heappop(heap)
        for pair in lo[1:]:
            pair[1] = '0' + pair[1]
        for pair in hi[1:]:
            pair[1] = '1' + pair[1]
        heapq.heappush(heap, [lo[0] + hi[0]] + lo[1:] + hi[1:])
    return sorted(heapq.heappop(heap)[1:], key=lambda p: (len(p[-1]), p))

# Função que auxilia na escrita dos bits
write_byte_buffer = 0
write_bits_in_buffer = 0
def write_bit(bit, flush=False):
    global write_byte_buffer, write_bits_in_buffer

    write_byte_buffer = (write_byte_buffer << 1) | bit
    write_bits_in_buffer += 1
    if write_bits_in_buffer == 8 or flush:
        value = write_byte_buffer
        write_bits_in_buffer = 0
        write_byte_buffer = 0
        return struct.pack('B', value)
    else:
        return None 

read_byte_buffer = 0
read_bits_in_buffer = 0
def read_bit(input):
    global read_byte_buffer, read_bits_in_buffer

    if read_bits_in_buffer == 0:
        read_byte_buffer = struct.unpack('B', input.read(1))[0]
        read_bits_in_buffer = 8
    bit = read_byte_buffer >> (read_bits_in_buffer - 1)
    read_bits_in_buffer -= 1
    return bit

# Função que escreve o cabeçalho do arquivo
def write_header(freqs, output):
    output.write(struct.pack('I', len(freqs)))
    output.write(struct.pack('I', img_width))
    output.write(struct.pack('I', img_height))
    output.write(struct.pack('?', is_gray))
    for sym, wt in freqs.items():
        #print(sym)
        output.write(struct.pack('B', sym))
        output.write(struct.pack('I', wt))

# Função que lê o cabeçalho do arquivo
def read_header(input):
    global img_width, img_height, is_gray

    freqs = {}
    num_syms = struct.unpack('I', input.read(4))[0]
    img_width = struct.unpack('I', input.read(4))[0]
    img_height = struct.unpack('I', input.read(4))[0]
    is_gray = struct.unpack('?', input.read(1))[0]
    for i in range(num_syms):
        sym = struct.unpack('B', input.read(1))[0]
        wt = struct.unpack('I', input.read(4))[0]
        freqs[sym] = wt
    return freqs

# Função que escreve os dados do arquivo
def write_data(array, tree, output):
    for sym in array:
        for bit in tree[sym][1]:
            byte = write_bit(int(bit))
            if byte is not None:
                output.write(byte)
    output.write(write_bit(0, flush=True))

# Função que lê os dados do arquivo
def read_data(input, tree):
    array = []
    codes = [code for sym, code in tree]

    for i in range(img_height * img_width * (1 if is_gray else 3)):
        byte = format(read_bit(input), 'b')
        while byte not in codes:
            byte += format(read_bit(input), 'b')
        array.append(byte)
    
    return array

# Função que comprime a imagem
def compress(array, output):
    freqs = calc_freq(array)
    tree = build_tree(freqs)
    write_header(freqs, output)
    write_data(array, tree, output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='huffman.py',
        description='Comprime e descomprime imagens usando o algoritmo de Huffman.',
        epilog='Exemplo de uso: python huffman.py -c -i imagem.png -o imagem.huff'
    )

    parser.add_argument('-c', '--compress', action='store_true', help='Comprime a imagem')
    parser.add_argument('-d', '--decompress', action='store_true', help='Descomprime a imagem')
    parser.add_argument('-i', '--input', type=str, help='Imagem de entrada', required=True)
    parser.add_argument('-o', '--output', type=str, help='Imagem de saída', required=True)
    parser.add_argument('-g', '--gray', action='store_true', help='Imagem em tons de cinza (padrão: RGB)')

    args = parser.parse_args()


    if args.compress:
        entrada = numpy.asarray(Image.open(args.input), numpy.uint8)

        shape = entrada.shape
        img_height = shape[0]
        img_width = shape[1]
        print(img_height, img_width)

        if args.gray:
            is_gray = True
            entrada = entrada.reshape((img_height * img_width))
        else:
            is_gray = False
            entrada = entrada.reshape((img_height * img_width * 3))

        with open(args.output, 'wb') as output:
            compress(entrada, output)

    elif args.decompress:
        with open(args.input, 'rb') as input:
            freqs = read_header(input)
            tree = build_tree(freqs)
            array = read_data(input, tree)
            if args.gray:
                array = numpy.asarray(array, numpy.uint8).reshape((img_height, img_width))
            else:
                array = numpy.asarray(array, numpy.uint8).reshape((img_height, img_width, 3))
 

        with open(args.output, 'wb') as output:
            Image.fromarray(array).save(output) 

    else:
        print('Nenhuma ação especificada. Use -c para comprimir ou -d para descomprimir.')