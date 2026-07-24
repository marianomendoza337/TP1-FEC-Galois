class GF:
    """
    Representa un campo de Galois GF(2^m).
    """
    def __init__(self, m, prim_poly):
        # El constructor recibe el orden m y el polinomio primitivo P(x)[cite: 1].
        self.m = m
        self.limit = 1 << m
        # prim_poly es un entero sin el término x^m (queda implícito)[cite: 1].
        self.prim_poly = prim_poly

    def add(self, a, b):
        """Suma de elementos del campo (operación XOR en binario)[cite: 1]."""
        return a ^ b

    def multiply(self, a, b):
        """Producto de elementos del campo con reducción por el polinomio primitivo[cite: 1]."""
        p = 0
        while b > 0:
            if b & 1:
                p ^= a
            b >>= 1
            a <<= 1
            # Reducción si el bit x^m se activa
            if a & self.limit:
                a ^= (self.limit | self.prim_poly)
        return p

    def power(self, a, n):
        """Potencia A^n para n > 0 usando exponenciación rápida[cite: 1]."""
        if n == 0:
            return 1
        res = 1
        base = a
        while n > 0:
            if n & 1:
                res = self.multiply(res, base)
            base = self.multiply(base, base)
            n >>= 1
        return res

    def inverse(self, a):
        """
        Inverso multiplicativo[cite: 1].
        Comportamiento para 0: Lanza una excepción (ZeroDivisionError)[cite: 1].
        En GF(2^m), a^(2^m - 2) es el inverso de a.
        """
        if a == 0:
            raise ZeroDivisionError("El inverso multiplicativo del elemento 0 no está definido.")
        return self.power(a, (1 << self.m) - 2)

    def divide(self, a, b):
        """
        División sobre el campo[cite: 1].
        Comportamiento para divisor 0: Lanza una excepción (ZeroDivisionError)[cite: 1].
        """
        if b == 0:
            raise ZeroDivisionError("División por el elemento 0 no definida.")
        return self.multiply(a, self.inverse(b))


class GFPoly:
    """
    Representa un polinomio con coeficientes en GF(2^m)[cite: 1].
    """
    def __init__(self, gf, coeffs):
        # Recibe una instancia de GF y una lista de coeficientes en orden decreciente de grado[cite: 1].
        self.gf = gf
        
        # Validar que cada coeficiente sea válido en el campo [0, 2^m - 1][cite: 1].
        for c in coeffs:
            if not (0 <= c < self.gf.limit):
                raise ValueError(f"El coeficiente {c} no es un elemento válido del campo GF(2^{self.gf.m}).")
                
        # Eliminar ceros a la izquierda, salvo en el polinomio nulo[cite: 1].
        while len(coeffs) > 1 and coeffs[0] == 0:
            coeffs = coeffs[1:]
        self.coeffs = list(coeffs) if coeffs else [0]

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __add__(self, other):
        """Suma de dos polinomios[cite: 1]."""
        max_len = max(len(self.coeffs), len(other.coeffs))
        c1 = [0] * (max_len - len(self.coeffs)) + self.coeffs
        c2 = [0] * (max_len - len(other.coeffs)) + other.coeffs
        res = [self.gf.add(a, b) for a, b in zip(c1, c2)]
        return GFPoly(self.gf, res)

    def __mul__(self, other):
        """Producto de dos polinomios[cite: 1]."""
        res = [0] * (len(self.coeffs) + len(other.coeffs) - 1)
        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(other.coeffs):
                prod = self.gf.multiply(a, b)
                res[i+j] = self.gf.add(res[i+j], prod)
        return GFPoly(self.gf, res)

    def __floordiv__(self, other):
        """División entera: retorna el cociente[cite: 1]."""
        q, _ = self._divmod(other)
        return q

    def __mod__(self, other):
        """División entera: retorna el resto[cite: 1]."""
        _, r = self._divmod(other)
        return r

    def _divmod(self, other):
        """Método interno para obtener cociente y resto[cite: 1]."""
        if len(other.coeffs) == 1 and other.coeffs[0] == 0:
            raise ZeroDivisionError("División de polinomios por el polinomio nulo.")
            
        dividend = self.coeffs[:]
        divisor = other.coeffs
        out_len = len(dividend) - len(divisor) + 1
        
        if out_len <= 0:
            return GFPoly(self.gf, [0]), GFPoly(self.gf, dividend)

        quotient = [0] * out_len
        for i in range(out_len):
            coef = self.gf.divide(dividend[i], divisor[0])
            quotient[i] = coef
            if coef != 0:
                for j in range(len(divisor)):
                    term = self.gf.multiply(coef, divisor[j])
                    dividend[i+j] = self.gf.add(dividend[i+j], term)
                    
        return GFPoly(self.gf, quotient), GFPoly(self.gf, dividend[out_len:])

    def scale(self, scalar):
        """Multiplicar todos los coeficientes por un escalar del campo[cite: 1]."""
        if not (0 <= scalar < self.gf.limit):
            raise ValueError("El escalar no pertenece al campo.")
        res = [self.gf.multiply(c, scalar) for c in self.coeffs]
        return GFPoly(self.gf, res)

    def evaluate(self, x):
        """Evaluación del polinomio en un punto x mediante el método de Horner[cite: 1]."""
        if not (0 <= x < self.gf.limit):
            raise ValueError("El punto a evaluar no pertenece al campo.")
        res = 0
        for c in self.coeffs:
            res = self.gf.add(self.gf.multiply(res, x), c)
        return res

    @classmethod
    def from_roots(cls, gf, roots):
        """Construcción a partir de un conjunto de raíces[cite: 1]."""
        poly = cls(gf, [1])
        for r in roots:
            # El término es (x - r). En GF(2^m), la resta es igual a la suma (XOR), por lo que es (x + r)
            root_poly = cls(gf, [1, r])
            poly = poly * root_poly
        return poly
