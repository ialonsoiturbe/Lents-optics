#include <stdio.h>
#include <math.h>

/*El programa ejecuta un bucle para cada rayo. En el bucle calcula la interseccion
del rayo con la primera lente, la intersección con la segunda lente y la interseccion 
con el eje (y=0). 
-La interseccion la he comprobado manualmente. a lapiz y papel. Funciona. Aunque solo
acepta valores positivos de h. Que se la va a hacer! Aunque no es inconveniente debido
a la simetria de la lente
-Para calcular la interseccion se utiliza el método de bisecciones. Mas comodo y con un 
error minimo.
-Para calcular la inclinacion que sufre el rayo se utilizan vectores. Se gira el vector del
normal a la superficie el angulo de refraccion. 
-Se imprime por pantalla y en un archivo datos.dat.
*/

/*/////////////////////////////////////DATOS////////////////////////////////////////////
R1 = +14.25	Radio de la primera lente
R2 = -9.47	Radio de la segunda lente
d  =  6.4	distancia entre lentes
n  = 1.61765	índice de refraccion
////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////CRITERIO DE SIGNOS/////////////////////////////////////
Como trabajamos con vectores es bastante necesario aclarar estas cosas:
ORIGEN DE REFERENCIA: Lo situamos en la primera lente. (Intersección lente-eje optico
SIGNOS PARA LOS VECTORES: Trabajamos en ejes cartesiano: 
				Hacia arriba 	 - positivo
				Hacia la derecha - positivo
SIGNOS PARA LOS ÁNGULOS: 
			giro antihorario 	- positivo
			giro horario		- negativo
Es bastante util, para implementar la Ley de Snell, tener este criterio de angulos
////////////////////////////////////////////////////////////////////////////////////////
*/
/*
==================================================================================
				MÉTODO DE BISECCIONES
==================================================================================
*/
double f(double t, double x0, double y0, double c, double R, double x) {
double f1;
f1 = pow(t*(x-x0) + y0, 2) + pow(x-c, 2) - R*R;
return f1;
}


double interseccion(double t,double x0,double y0,double c,double R) {

	int N;
	double x1, x2;
	double f1, f3;
	double x3, f13;

	// limites del intervalo inicial. Deben estar en los límites de la lente.
	x1 = 0.0;
	x2 = 6.4;
	// Inicializamos a 0 el contador de iteraciones
	N=0;
	

	// Bucle do while del metodo de biseccion 
	do {

	// Calculamos punto medio del intervalo
	// y los valores f1, f2, f3
	x3 = 0.5*(x2+x1);
	f1 = f(t, x0, y0, c, R, x1);
	f3 = f(t, x0, y0, c, R, x3);

	// Dependiendo de los signos de f(v1) y f(v3), elegimos 
	// los limites del nuevo intervalo  
	f13 = f1*f3;
	if(f13<0) {
	x1 = x1;
	x2 = x3;} 
	else {
	x1 = x3; 
	x2 = x2;
	}

	// Incrementamos el contador de iteraciones en 1
	N = N+1;

		// Paramos el bucle para mas de 1000 iteraciones
		if (N>1000) {
		printf("Numero maximo de Iteraciones alcanzado \n");
		return 0;
		}

	// El bucle se para si el valor absoluto de f(V) es menor que 0.000001
	} while (fabs(f3)>0.000001);
	
	return x3;

	}
/*
==================================================================================
				OTRAS FUNCIONES 
==================================================================================
*/
void giro (double U[], double V[], double tita) //(Ux, Uy) vector a girar un ángulo TITA en sentido antihorario.
	{
	V[0] = U[0]*cos(tita) - U[1]*sin(tita);
	V[1] = U[0]*sin(tita) + U[1]*cos(tita);
	}

double angulo_vectores(double V[], double U[]) {
	double modulo, angulo;
	modulo = sqrt(V[0]*V[0] + V[1]*V[1]);
	V[0] = V[0]/modulo; 
	V[1] = V[1]/modulo;

	modulo = sqrt(U[0]*U[0] + U[1]*U[1]);
	U[0] = U[0]/modulo; 
	U[1] = U[1]/modulo;

	angulo = acos ( V[0] * U[0] + V[1]*U[1]);
return angulo;
}

//================================================================================
// 					PROGRAMA 
//================================================================================

int main () {
	double hmax, hpaso, h, N;
	long i;
	long j;

	///////////////////Variables para el bucle de rayos///////////////////////
	double R1 = 14.25;
	double R2 = 9.47; //Radios
	double d  = 6.4; //Distancia de separación de lentes
	double n = 1.61765; //Indice de refracción
	
	double p1[2], p2[2], p3[2];
	double xaux, yaux, t, m;
	double V[2], U[2];
	double alpha1, alpha2, ep1, ep2;

	//////////////////////////////////////////////////////////////////////////

	printf("Indique el radio (mm) de apertura dentro del intervalo (0, 6.3)\n");
	scanf("%lf", &hmax);
	printf("Número de rayos a trazar:\n");
	scanf("%lf", &N); 

	hpaso = hmax / N;

	FILE *output;
	output = fopen("datos.dat", "w");
	/////Cálculo de rayos///////
	/*Contaremos el número de rayos (i). La altura de cada rayo será
	h = i*hpaso. No calculamos el h = 0 por ser trivial. El valor máximo 
	será N*hpaso = hmax;
	*/
	
	for (i = 1; i<= (long)N; i++)
	{
		printf("-------%ld-----------\n", i);
		h = i*hpaso; 
		for (j= 0; j<2; j++) //"vaciamos" los arrays en cada iteración. Mejor que no den problemas
		{
			p1[j] = 0;
			p2[j] = 0;
			p3[j] = 0;
			V[j] = 0;
			U[j] =0;
		}
		
		//Puntos:
		//1)
		xaux = interseccion(0, -1, h, R1, R1);
	 	yaux = h;
		p1[0] = xaux;
		p1[1] = yaux;

		//2)
		V[0] = 1.0;	//Para el criterio de signos. El ángulo formado para aplicar Snell debe ser 
		V[1] = 0;	//el comprendido entre este vector y el radio	     
		U[0] = R1 - p1[0];//Vector radio. Considerado así, nos será fácil girarlo de forma correcta
		U[1] = -1.0*p1[1];

		ep1 = angulo_vectores( V, U); //Nuestro ángulo tiene que ser positivo según los signos de antes
					      //Aquí no hay limitación de ángulo. je.
		ep2 = asin( sin(ep1) / n ); //Ley de Snell (ÁNGULO: SALE POSITIVO)
		giro ( U, V, ep2); //Gira el vector radio en un angulo ep2 
		t = V[1]/V[0]; //Pendiente: t = dy/dx;
		xaux = interseccion(t, p1[0], p1[1], d-R2, R2);
		yaux = t*(xaux - p1[0]) + p1[1];
		p2[0] = xaux;
		p2[1] = yaux;

		//3)
		U[0] = p2[0] - (d- R2); //Vector radio
		U[1] = p2[1]; 
		//El vector V está definido de antes
		alpha1 = angulo_vectores( V, U);//Angulo entre radio y recta 
		alpha1 = -1.0*alpha1;	//El ángulo es negativo.				  
		if ( fabs(alpha1) > 0.666426058) {
			printf("Ángulo límite superado: reflexión total. aplha1 >0.666426058 \n");
			printf("p1 (%lf, %lf)", p1[0], p1[1]);
			printf("t = %lf, ep1 = %lf, ep2 = %lf\n", t, ep1, ep2);
			fprintf(output,"%lf		%.8lf\n", h, p3[0]);			
			fclose(output);
			return 0;
			}
		alpha2 = asin( sin(alpha1) * n );	//Ley de Snell (NEGATIVO)
		giro ( U, V, alpha2);
		m = V[1]/V[0];
		p3[0] = (1.0/m) * (0 - p2[1]) + p2[0]; //x = 1/m *(y-y0) + x0
		p3[1] = 0;
		fprintf(output,"%lf		%.8lf\n", h, p3[0]);
	}

	printf("FIN DEL PROGRAMA\n");
	fclose(output);
	return 0;
}
