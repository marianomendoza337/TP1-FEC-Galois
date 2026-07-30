from typing import List, Tuple, Dict

class GaloisField:
    """
    Representa un campo de Galois GF(2^m).
    Los elementos se representan como enteros en el rango [0, 2^m - 1].
    """
    def __init__(self, m: int, prim_poly: int):
        """
        Inicializa el campo de Galois.

        Args:
            m: Orden del campo (m bits).
            prim_poly: Polinomio primitivo representado como un entero de m bits
                       (sin el término x^m, que queda implícito).
        """
        self.m = m
        self.prim_poly = prim_poly
        self.order = 1 << m   #Calcula el numero total de elementos = 2^m

    #Definimos los metodos de la clase

    def _check_valid(self, a: int) -> None:
        """Verifica que el elemento pertenezca al campo."""
        if not (0 <= a < self.order):
            raise ValueError(f"El elemento {a} no pertenece al rango válido [0, {self.order - 1}].")   #Lanza un error si no esta comprendido en el rango del campo

    def add(self, a: int, b: int) -> int:
        """Suma de dos elementos en el campo (XOR bit a bit)."""
        self._check_valid(a)
        self._check_valid(b)
        return a ^ b    #En py el operador ^ representa la funcion XOR bit a bit

    def mul(self, a: int, b: int) -> int:
        """Multiplicación de dos elementos en el campo con reducción polinomial."""
        self._check_valid(a)
        self._check_valid(b)
        res = 0
        for _ in range(self.m):
            if b & 1:     #el bit de la derecha es 1?
                res = res ^ a  # Sumamos a usando XOR
            # Detectamos si el bit más significativo (x^(m-1)) está encendido
            carry = a & (1 << (self.m - 1))           #operador << desplazamiento a la izquierda
            # Desplazamos y aplicamos máscara para mantener m bits
            a = (a << 1) & (self.order - 1)
            # Si hubo desbordamiento, aplicamos reducción con el polinomio primitivo
            if carry:
                a ^= self.prim_poly
            b >>= 1
        return res

    def power(self, a: int, n: int) -> int:
        """
        Potencia A^n para n > 0 utilizando exponenciación rápida.
        Al usar esta forma de exponenciacion rapida, supongamos 2**8=256
        De la forma tradicional tardariamos 256 iteraciones, en cambio con este
        metodo solo 8 iteraciones.

        """
        self._check_valid(a)
        if n < 0:
            raise ValueError("El exponente debe ser n >= 0.")
        if n == 0:
            return 1

        res = 1 #variable donde se guardaran los numeros
        base = a    # ahora tenemos a**1
        while n > 0:    #mientras n sea positivo
            if n & 1:   #a nivel bits, el numero n es impar?
                res = self.mul(res, base) #multiplica el valor que hay en res por base = a
            base = self.mul(base, base)   #Elevamos la base al cuadrado
            n >>= 1        #Divide el exponente por 2 e ignora los decimales, ejemplo 5/2=2
        return res          #Devuelve el resto

    def inv(self, a: int) -> int:
        """
        Inverso multiplicativo de un elemento a.
        Lanza ZeroDivisionError si el elemento es 0.
        """
        self._check_valid(a)
        if a == 0:
            raise ZeroDivisionError("El elemento 0 no tiene inverso multiplicativo.")
        # En GF(2^m), a^(2^m - 1) = 1, por ende a^(-1) = a^(2^m - 2)
        return self.power(a, self.order - 2)    #Esto sale del teorema de Fermat Si m=4, entonces hay 2**4=16, por Fermat  a^(2^m - 1) = 1

    def div(self, a: int, b: int) -> int:
        """
        División de a entre b.
        Lanza ZeroDivisionError si el divisor b es 0.
        """
        self._check_valid(a)
        self._check_valid(b)
        if b == 0:
            raise ZeroDivisionError("División por cero en el Campo de Galois.")
        return self.mul(a, self.inv(b))     #Basicamente en lugar de dividir por b, multiplicamos a por el inverso multiplicativo de b

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GaloisField):  #Si el otro objeto no pertenece al campo de Galois?
            return False
        return self.m == other.m and self.prim_poly == other.prim_poly    #Los objetos son iguales


class GFPoly:
    """
    Representa un polinomio con coeficientes en GF(2^m).
    """
    def __init__(self, gf: GaloisField, coeffs: List[int]):
        """
        Inicializa el polinomio.

        Args:
            gf: Instancia de la clase GaloisField.
            coeffs: Lista de coeficientes en orden decreciente de grado.
        """
        self.gf = gf

        # Validar que cada coeficiente pertenezca al campo
        for c in coeffs:
            self.gf._check_valid(c)     #Evaluo que todos los coeficientes pertenescan al campo

        # Remover ceros a la izquierda (excepto para el polinomio nulo)
        while len(coeffs) > 1 and coeffs[0] == 0: #Mientras la cantidad de coeficientes sea mayor a 1 y ademas el primer numero de la lista sea un 0
            coeffs = coeffs[1:]   #Slicing: Reemplaza la lista actual por una copia que empiece desde la posición 1 en adelante

        if not coeffs:  #Si la lista esta vacia
            coeffs = [0]  

        self.coeffs = coeffs  #Guardo en memoria

    @property       #Con decorador property evitamos escribir la funcion con ()
    def degree(self) -> int:  
        """Devuelve el grado del polinomio. El polinomio nulo tiene grado -1."""
        if len(self.coeffs) == 1 and self.coeffs[0] == 0:   #¿La lista tiene un solo elemento y ese solo elemento es cero?
            return -1   #Devuelvo -1
        return len(self.coeffs) - 1   #Como en la funcion anterior la lista ya esta limpia, la cantidad de elementos -1 nos indica el grado del polinomio

    def __add__(self, other: 'GFPoly') -> 'GFPoly':
        """Suma de dos polinomios."""
        max_len = max(len(self.coeffs), len(other.coeffs))  #Indico cual es la lista mas larga
        c1 = [0] * (max_len - len(self.coeffs)) + self.coeffs #Agrego una lista con ceros al polinomio mas pequeño
        c2 = [0] * (max_len - len(other.coeffs)) + other.coeffs #Si el polinomio es el mas grande sumo una lista vacia, basicamente no pasa nada
        new_coeffs = [self.gf.add(a, b) for a, b in zip(c1, c2)]      #Con zip junto los coeficientes del mismo grado en pareas, recorro con el for esos pares y luego uso XOR para sumar
        return GFPoly(self.gf, new_coeffs)  #Nueva lista con la suma

    def __mul__(self, other: 'GFPoly') -> 'GFPoly':   
        """Producto de dos polinomios."""
        if self.degree == -1 or other.degree == -1:     #Si alguno de los polinomios es el polinomio nulo
            return GFPoly(self.gf, [0])     #Devuelvo el polinomio nulo

        res_degree = self.degree + other.degree   #Sumo el grado de los polinomios para saber el grado final
        res_coeffs = [0] * (res_degree + 1)       #Creo una lista de de ceros de tamñano del grado +1

        for i, a in enumerate(self.coeffs):       #i, j son las posiciones de las listas y a,b los valores de las listas
            for j, b in enumerate(other.coeffs):
                prod = self.gf.mul(a, b)          #Usamos la multiplicacion antes definida
                res_coeffs[i + j] = self.gf.add(res_coeffs[i + j], prod)    #Ve a la posicion i+j, usa self.gf.add para sumar lo que ya habia

        return GFPoly(self.gf, res_coeffs)        #Devuelve el polinomio resultante

    def __divmod__(self, other: 'GFPoly') -> Tuple['GFPoly', 'GFPoly']:
        """
        División entera de dos polinomios.
        Retorna (cociente, resto).
        """
        if other.degree == -1:
            raise ZeroDivisionError("División polinomial por el polinomio nulo.") #Devuelve error si el divisor es el polinomio nulo

        if self.degree < other.degree:      
            return GFPoly(self.gf, [0]), GFPoly(self.gf, self.coeffs)   #El cociente es cero y el resto basicamente es el dividendo original

        rem_coeffs = self.coeffs.copy()           #Coeficientes del resto, al inicio tiene la misma cantidad de elementos que el dividendo
        out_degree = self.degree - other.degree   #Resta de exponentes
        quot_coeffs = [0] * (out_degree + 1)      #Creo una lista con ceros para ir completando el cociente

        inv_lead = self.gf.inv(other.coeffs[0])   #Calculo el inverso multiplicativo del mayor coeficiente del divisor

        for i in range(out_degree + 1):
            if rem_coeffs[i] != 0:                #Si el numero no es cero se procede a eliminarlo
                mult = self.gf.mul(rem_coeffs[i], inv_lead) 
                quot_coeffs[i] = mult

                for j in range(len(other.coeffs)):
                    prod = self.gf.mul(mult, other.coeffs[j])
                    rem_coeffs[i + j] = self.gf.add(rem_coeffs[i + j], prod)

        return GFPoly(self.gf, quot_coeffs), GFPoly(self.gf, rem_coeffs)

    def __floordiv__(self, other: 'GFPoly') -> 'GFPoly':
        """Cociente de la división entera de polinomios."""
        q, _ = self.__divmod__(other)     #Aplico la division, pero solo me quedo con el cociente, no importa el resto por eso el _
        return q                           #Devuelvo solo el cociente de la division entera

    def __mod__(self, other: 'GFPoly') -> 'GFPoly':
        """Resto de la división entera de polinomios."""
        _, r = self.__divmod__(other)           
        return r                            #Devuelve el resto de la division entera

    def scale(self, scalar: int) -> 'GFPoly':
        """Multiplica todos los coeficientes del polinomio por un escalar del campo."""
        self.gf._check_valid(scalar)          #Verificamos que el escalar pertenesca al campo de Galois
        new_coeffs = [self.gf.mul(c, scalar) for c in self.coeffs]      #Recorro todos los coeficientes de coeffs, y multiplico por el escalar
        return GFPoly(self.gf, new_coeffs)      #Con la lista nueva de coeficientes crea un objeto

    def evaluate(self, x: int) -> int:
        """
        Evaluación del polinomio en un punto x utilizando el método de Horner.
        """
        self.gf._check_valid(x)     #verificamos que x pertenesca al campo de Galois
        res = 0                     #Resultado con valor inicial 0
        for c in self.coeffs:       #Recorremos los coeficientes
            res = self.gf.add(self.gf.mul(res, x), c)     # res * x + c   -> Metodo de Horner
        return res

    @classmethod
    def from_roots(cls, gf: GaloisField, roots: List[int]) -> 'GFPoly':   #cls es clase y es para no escribir GFPoly todo el tiempo
        #roots es la lista de raices
        """
        Construye el polinomio dado un conjunto de raíces.
        Como en GF(2^m) sumar y restar es lo mismo (-r = r), (x - r) es (x + r).
        """
        poly = cls(gf, [1])   #Se crea un polinomio unitario
        for r in roots:       
            gf._check_valid(r) #vericar que la raiz se encuentre dentro del campo
            term = cls(gf, [1, r])  #Crea un polinomio usando el 1 y la raiz
            poly = poly * term      
        return poly                 #Devuelve el polinomio

    def __eq__(self, other: object) -> bool:      #Definimos que __eq__ es para comparar polinomios "=="
        if not isinstance(other, GFPoly):         #El otro objeto es un polinomio
            return False
        return self.coeffs == other.coeffs and self.gf == other.gf      #Los polinomios son iguales si poseen los mismos coeficientes y viven en el mismo campo

    def __repr__(self) -> str:
        return f"GFPoly({self.coeffs})"

# ============================================================================
# FUNCIONES AUXILIARES PARA REPRESENTACIÓN
# ============================================================================

def build_alpha_map(gf: GaloisField, alpha: int = 2) -> Tuple[Dict[int, str], Dict[str, int]]:
    """Crea un mapa bidireccional entre enteros y potencias de alpha."""
    a_map = {0: "0"}      #Diccionario de numeros a letras
    val_map = {"0": 0}    #Diccionario de letras a numeros
    for i in range(gf.order - 1):
        val = gf.power(alpha, i)
        label = f"a^{i}" if i > 1 else ("a" if i == 1 else "1")
        a_map[val] = label
        val_map[label] = val
    return a_map, val_map

def to_a(val: int, a_map: Dict[int, str]) -> str:
    """Convierte un valor entero a su representación en potencia de alpha."""
    return a_map.get(val, "?")
