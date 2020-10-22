//Headers del programa
#ifndef VSMAT_H
#define VSMAT_H

/*
 * Creamos una estructura para un vector de 2 coordenadas.
 */
typedef struct {
	/* data */
	float w;
	float x;
} vec2; 

/*
 * Creamos una estructura para un vector de 3 coordenadas.
 */
typedef struct {
	/* data */
	float w;
	float x;
	float y;
} vec3;

/*
 * Creamos una estructura para un vector de 4 coordenadas.
 */
typedef struct {
	/* data */
	float w;
	float x;
	float y;
	float z;
} vec4;

/*
 * Creamos una estructura para una matriz de 4x4 a partir de un arreglo unidimensional.
 */
typedef float mat4[16];
mat4 m;

#endif