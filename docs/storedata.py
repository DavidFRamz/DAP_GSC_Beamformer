import pickle
def save(filename, *args, glob):
    # Obtener el diccionario de var globales
    #glob = globals()
    d = {}
    for v in args:
        # Copiar los valores al diccionario
        d[v] = glob[v]
    with open(filename, 'wb') as f:
        # Ponerlos en un archivo
        pickle.dump(d, f)

def load(filename):
    # Obtener el diccionario de vars globales
    glob = globals()
    with open(filename, 'rb') as f:
        for k, v in pickle.load(f).items():
            # Definir y Asiganr c/var global 
            glob[k] = v
