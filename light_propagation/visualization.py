import numpy as np

def main():
    with open("binary.data", 'rb') as f:
        rectype = np.dtype(np.float64)
        data = np.fromfile(f, dtype=rectype)

    np.savetxt("test.txt", data)

if __name__ == '__main__':
    main()
