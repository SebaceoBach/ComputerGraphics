/*
 *
 * Implementación de algoritmos para el manejo de vectores de 2,3 y 4 coordenadas y
 * matrices de 4x4 con arreglos unidimensionales
 * @author Sebastian Salazar
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vsmat.h"

/*
 * Reliza el producto punto entre dos vectores u,v y devuelve un flotante.
 */
float vsDot(const vec3 *v, const vec3 *u) {
	float s = v->w*u->w;
	float t = v->x*u->x;
	float r = v->y*u->y;

	return (s+t+r);
}
/*
 * Realiza el producto cruz entre dos vectores u,v y guarda el resultado en dest.
 */
void vsCross(vec3 *dest, const vec3 *v, const vec3 *u){
	float i = (v->x*u->y)-(v->y*u->x);
	float j = (v->w*u->y)-(v->y*u->w);
	float k = (v->w*u->x)-(v->x*u->w);

	dest -> w = i;
	dest -> x = j;
	dest -> y = k;	
}
/*
 * Normaliza un vector
 */
void vsNormalize(vec3 *v){
	sqrtf((v->w*v->w) + (v->x*v->x) + (v->y*v->y));
}

/*
 * Hace la identidad de una matriz
 */
void vsIdentity(mat4 m){
	
	m[0] = 1;
	m[1] = 0;
	m[2] = 0;
	m[3] = 0;
	m[4] = 0;
	m[5] = 1;
	m[6] = 0;
	m[7] = 0;
	m[8] = 0;
	m[9] = 0;
	m[10] = 1;
	m[11] = 0;
	m[12] = 0;
	m[13] = 0; 
	m[14] = 0;
	m[15] = 1;
}


/*
 * Realiza el producto entre dos matrices m1*m2 y guarda el resultado en dest
 */
void vsMatMul(mat4 dest, const mat4 m1, const mat4 m2){
	dest[0] = (m1[0]*m2[0]) + (m1[4]*m2[1]) + (m1[8]*m2[2]) + (m1[12]*m2[3]);
	dest[4] = (m1[0]*m2[4]) + (m1[4]*m2[5]) + (m1[8]*m2[6]) + (m1[12]*m2[7]);
	dest[8] = (m1[0]*m2[8]) + (m1[4]*m2[9]) + (m1[8]*m2[10]) + (m1[12]*m2[11]);
	dest[12] = (m1[0]*m2[12]) + (m1[4]*m2[13]) + (m1[8]*m2[14]) + (m1[12]*m2[15]);
	//Cambio de fila
	dest[1] = (m1[1]*m2[0]) + (m1[5]*m2[1]) + (m1[9]*m2[2]) + (m1[13]*m2[3]);
	dest[5] = (m1[1]*m2[4]) + (m1[5]*m2[5]) + (m1[9]*m2[6]) + (m1[13]*m2[7]);
	dest[9] = (m1[1]*m2[8]) + (m1[5]*m2[9]) + (m1[9]*m2[10]) + (m1[13]*m2[11]);
	dest[13] = (m1[1]*m2[12]) + (m1[5]*m2[13]) + (m1[9]*m2[14]) + (m1[13]*m2[15]);
	//Cambio de fila
	dest[2] = (m1[2]*m2[0]) + (m1[6]*m2[1]) + (m1[10]*m2[2]) + (m1[14]*m2[3]);
	dest[6] = (m1[2]*m2[4]) + (m1[6]*m2[5]) + (m1[10]*m2[6]) + (m1[14]*m2[7]);
	dest[10] = (m1[2]*m2[8]) + (m1[6]*m2[9]) + (m1[10]*m2[10]) + (m1[14]*m2[11]);
	dest[14] = (m1[2]*m2[12]) + (m1[6]*m2[13]) + (m1[10]*m2[14]) + (m1[14]*m2[15]);
	//Cambio de fila
	dest[3] = (m1[3]*m2[0]) + (m1[7]*m2[1]) + (m1[11]*m2[2]) + (m1[15]*m2[3]);
	dest[7] = (m1[3]*m2[4]) + (m1[7]*m2[5]) + (m1[11]*m2[6]) + (m1[15]*m2[7]); 
	dest[11] = (m1[3]*m2[8]) + (m1[7]*m2[9]) + (m1[11]*m2[10]) + (m1[15]*m2[11]);
	dest[15] = (m1[3]*m2[12]) + (m1[7]*m2[13]) + (m1[11]*m2[14]) + (m1[15]*m2[15]);
}
/*
 * Realiza la transpuesta de una matriz m y guarda la matriz resultante en dest.
 */
void vsTranspose(mat4 dest, mat4 m){
	dest[0] = m[0];
	dest[1] = m[4];
	dest[2] = m[8];
	dest[3] = m[12];
	//Cambio de columna
	dest[4] = m[1];
	dest[5] = m[5];
	dest[6] = m[9];
	dest[7] = m[13];
	//Cambio de columna
	dest[8] = m[2];
	dest[9] = m[6];
	dest[10] = m[10];
	dest[11] = m[14];
	//Cambio de columna
	dest[12] = m[3];
	dest[13] = m[7];
	dest[14] = m[11];
	dest[15] = m[15];
}

/*
 * Realiza la matriz de traslación de m.
 */
void vsTranslate(mat4 m, float x, float y, float z){
	m[0] = (m[0]*1) + (m[4]*0) + (m[8]*0) + (m[12]*0);
	m[4] = (m[0]*0) + (m[4]*1) + (m[8]*0) + (m[12]*0);
	m[8] = (m[0]*0) + (m[4]*0) + (m[8]*1) + (m[12]*0);
	m[12] = (m[0]*x) + (m[4]*y) + (m[8]*z) + (m[12]*1);
	//Cambio de fila
	m[1] = (m[1]*1) + (m[5]*0) + (m[9]*0) + (m[13]*0);
	m[5] = (m[1]*0) + (m[5]*1) + (m[9]*0) + (m[13]*0);
	m[9] = (m[1]*0) + (m[5]*0) + (m[9]*1) + (m[13]*0);
	m[13] = (m[1]*x) + (m[5]*y) + (m[9]*z) + (m[13]*1);
	//Cambio de fila
	m[2] = (m[2]*1) + (m[6]*0) + (m[10]*0) + (m[14]*0);
	m[6] = (m[2]*0) + (m[6]*1) + (m[10]*0) + (m[14]*0);
	m[10] = (m[2]*0) + (m[6]*0) + (m[10]*1) + (m[14]*0);
	m[14] = (m[2]*x) + (m[6]*y) + (m[10]*z) + (m[14]*1);
	//Cambio de fila
	m[3] = (m[3]*1) + (m[7]*0) + (m[11]*0) + (m[15]*0);
	m[7] = (m[3]*0) + (m[7]*1) + (m[11]*0) + (m[15]*0); 
	m[11] = (m[3]*0) + (m[7]*0) + (m[11]*1) + (m[15]*0);
	m[15] = (m[3]*x) + (m[7]*y) + (m[11]*z) + (m[15]*1);
}

/*
 * Realiza la matriz de rotación de m.
 */
void vsRotate(mat4 m, float angle, float x, float y, float z){
	float s = sinf(angle);
	float c = cosf(angle);

	/*
	* Definimos la matriz R
	*/
	mat4 r = {(((x*x)*(1-c))+c),
		(((y*x)*(1-c))+(z*s)),
		(((z*x)*(1-c))-(y*s)),
		0,
		//Cambio de columna
		(((x*y)*(1-c))-(z*s)),
		(((y*y)*(1-c))+c),
		(((z*y)*(1-c))+(x*s)),
		0,
		//Cambio de columna
		(((x*z)*(1-c))+(y*s)),
		(((y*z)*(1-c))-(x*s)),
		(((z*z)*(1-c))+c),
		0,
		//Cambio de columna
		0,0,0,1};

		//Multiplicamos m*r y guardamos el resultado en la misma matriz m
		vsMatMul(m,m,r); 

}

/*
 * Reliza el producto entre m y la matriz de homotecia.
 */
void vsScale(mat4 m, float x, float y, float z) {
	mat4 s = {x,0,0,0,
			  0,y,0,0,
			  0,0,z,0,
			  0,0,0,1};
	vsMatMul(m,m,s);

}
/*
 * Reliza el producto entre m y la matriz de proyección ortogonal.
 */
void vsOrtho(mat4 m, float left, float right, float bottom, 
	float top, float near,float far) {

	float tx = -((right+left)/(right-left));
	float ty = -((top+bottom)/(top-bottom));
	float tz = -((far+near)/(far-near));

	mat4 o = {2/(right-left),0,0,0,
			  0, 2/(top-bottom),0,0,
			  0,0, (-2)/(far-near),0,
			  tx,ty,tz,1};

	vsMatMul(m,m,o);

}

/*
 * Realiza el producto entre m y la matriz de proyección en perspectiva.
 */
void vsPerspective(mat4 m, float fovy, float aspect, float znear, float zfar) {

	float f = tanf(fovy/2); //revisar si es cotanf o tanf

	mat4 p = {f/aspect,0,0,0,
			  0,f,0,0,
			  0,0, (zfar+znear)/(znear-zfar), -1,
			  0,0, (2*zfar*znear)/(znear-zfar),0};

	vsMatMul(m,m,p);
	
}

/*
 * Creamos una matriz v a partir del vector unitario u
 * del producto cruz entre s y f
 */
void vsLookAt(mat4 m, float eyex, float eyey, float eyez, float centx, float centy,
float centz, float upx, float upy, float upz) {
	
	vec3 *f;
	vec3 *up;
	vec3 *s;
	vec3 *u;

	f -> w = (centx-eyex);
	f -> x = (centy-eyey);
	f -> y = (centz-eyez);

	up -> w = upx; 
	up -> x = upy;
	up -> y = upz;

	vsCross(s,f,up);
	vsNormalize(s);
	vsCross(u,s,f);

	mat4 v = {s -> w, u -> w, -(f -> w), 0,
			  s -> x, u -> x, -(f -> x), 0,
			  s -> y, u -> y, -(f -> y), 0,
			  0,0,0,1};

	vsMatMul(m,m,v);
	vsTranslate(m, -eyex, -eyey, -eyez);
}

/*
 * Realiza el determinante de la matriz m.
 */
float vsDeterminant(mat4 m){
	return (m[3] * m[6]  * m[9]  * m[12]   -  m[2] * m[7] * m[9]  * m[12]   -
         m[3] * m[5]  * m[10] * m[12]   +  m[1] * m[7] * m[10] * m[12]   +
         m[2]  * m[5]  * m[11] * m[12]   -  m[1] * m[6]  * m[11] * m[12]   -
         m[3] * m[6]  * m[8]  * m[13]   +  m[2] * m[7] * m[8]  * m[13]   +
         m[3] * m[4]  * m[10] * m[13]   -  m[0] * m[7] * m[10] * m[13]   -
         m[2]  * m[4]  * m[11] * m[13]   +  m[0] * m[6]  * m[11] * m[13]   +
         m[3] * m[5]  * m[8]  * m[14]  -  m[1] * m[7] * m[8]  * m[14]  -
         m[3] * m[4]  * m[9]  * m[14]  +  m[0] * m[7] * m[9]  * m[14]  +
         m[1]  * m[4]  * m[11] * m[14]  -  m[0] * m[5]  * m[11] * m[14]  -
         m[2]  * m[5]  * m[8]  * m[15]  +  m[1] * m[6]  * m[8]  * m[15]  +
         m[2]  * m[4]  * m[9]  * m[15]  -  m[0] * m[6]  * m[9]  * m[15]  -
         m[1]  * m[4]  * m[10] * m[15]  +  m[0] * m[5]  * m[10] * m[15]);
} 

/*
 * Realiza la matriz inversa de m.
 */
void vsInverse(mat4 dest, mat4 m) {
	double inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
       

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        dest[i] = inv[i] * det;

    


}

int main(int argc, char const *argv[]) {
	/* code */

	printf("\n");
	return 0;
}