# Galois Field & Polynomial Arithmetic en Python 

Una implementación orientada a objetos, pura y eficiente de la aritmética de Campos de Galois $GF(2^m)$ y operaciones con polinomios. 

Este proyecto proporciona un "motor matemático" que reemplaza la aritmética tradicional por matemáticas de campos finitos, sentando las bases algorítmicas para sistemas de corrección de errores.

##  Características Principales

*   **Aritmética de Campo Finito $GF(2^m)$:** Suma (XOR), multiplicación con reducción polinomial, y cálculo de inversos multiplicativos.
*   **Exponenciación Rápida:** Implementación eficiente de potencias en tiempo logarítmico.
*   **Polimorfismo y Sobrecarga de Operadores:** Operaciones matemáticas naturales con polinomios usando símbolos de Python (`+`, `*`, `//`, `%`, `==`).
*   **División Larga de Polinomios:** Soporte completo para calcular el cociente y el resto (`divmod`).
*   **Evaluación Optimizada:** Uso del **Método de Horner** para evaluar polinomios de manera eficiente sin calcular potencias explícitas.
*   **Generación desde Raíces:** Constructor especializado (`@classmethod`) para ensamblar polinomios complejos a partir de sus raíces.

##  Requisitos e Instalación

Este proyecto está escrito en Python y no tiene dependencias externas. 
Requiere **Python 3.6+** (debido al uso de Type Hints y f-strings).

