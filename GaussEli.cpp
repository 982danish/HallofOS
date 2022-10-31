// C# program to demonstrate working of Gaussian Elimination method
using System;

class GFG{
	
public static int N = 3;

static void gaussianElimination(double [,]mat)
{

	/* reduction into r.e.f. */
	int singular_flag = forwardElim(mat);
	
	/* if matrix is singular */
	if (singular_flag != -1)
	{
		Console.WriteLine("Singular Matrix.");
		
		/* if the RHS of equation corresponding to
		zero row is 0, * system has infinitely
		many solutions, else inconsistent*/
		if (mat[singular_flag,N] != 0)
			Console.Write("Inconsistent System.");
		else
			Console.Write("May have infinitely " +
						"many solutions.");
		
		return;
	}
	
	/* get solution to system and print it using
	backward substitution */
	backSub(mat);
}

static void swap_row(double [,]mat, int i, int j)
{
	
	
	for(int k = 0; k <= N; k++)
	{
		double temp = mat[i, k];
		mat[i, k] = mat[j, k];
		mat[j, k] = temp;
	}
}
	
static void print(double [,]mat)
{
	for(int i = 0; i < N; i++, Console.WriteLine())
		for(int j = 0; j <= N; j++)
			Console.Write(mat[i, j]);
			
	Console.WriteLine();
}

static int forwardElim(double [,]mat)
{
	for(int k = 0; k < N; k++)
	{
		
		int i_max = k;
		int v_max = (int)mat[i_max, k];
		
		/* find greater amplitude for pivot if any */
		for(int i = k + 1; i < N; i++)
		{
			if (Math.Abs(mat[i, k]) > v_max)
			{
				v_max = (int)mat[i, k];
				i_max = i;
			}
		
			/* If a principal diagonal element is zero,
			* it denotes that matrix is singular, and
			* will lead to a division-by-zero later. */
			if (mat[k, i_max] == 0)
				return k; // Matrix is singular
			
			/* Swap the greatest value row with
			current row
			*/
			if (i_max != k)
				swap_row(mat, k, i_max);
			
			for(int i = k + 1; i < N; i++)
			{
			
				/* factor f to set current row kth element
				* to 0, and subsequently remaining kth
				* column to 0 */
				double f = mat[i, k] / mat[k, k];
				
				/* subtract fth multiple of corresponding
				kth row element*/
				for(int j = k + 1; j <= N; j++)
					mat[i, j] -= mat[k, j] * f;
				
				/* filling lower triangular matrix with
				* zeros*/
				mat[i, k] = 0;
			}
		}
	}
	
	return -1;
}

static void backSub(double [,]mat)
{
	
	double []x = new double[N];
	
	/* Start calculating from last equation up to the
	first */
	for(int i = N - 1; i >= 0; i--)
	{
	
		/* start with the RHS of the equation */
		x[i] = mat[i,N];
		
		/* Initialize j to i+1 since matrix is upper
		triangular*/
		for(int j = i + 1; j < N; j++)
		{
		
			/* subtract all the lhs values
			* except the coefficient of the variable
			* whose value is being calculated */
			x[i] -= mat[i,j] * x[j];
		}
		
		/* divide the RHS by the coefficient of the
		unknown being calculated */
		x[i] = x[i] / mat[i,i];
	}
	
	Console.WriteLine();
	Console.WriteLine("Solution for the system:");
	for(int i = 0; i < N; i++)
	{
		Console.Write("{0:F6}", x[i]);
		Console.WriteLine();
	}
}

public static void Main(String[] args)
{
	
	/* input matrix */
	double [,]mat = { { 3.0, 2.0, -4.0, 3.0 },
					{ 2.0, 3.0, 3.0, 15.0 },
					{ 5.0, -3, 1.0, 14.0 } };
	
	gaussianElimination(mat);
}
}

