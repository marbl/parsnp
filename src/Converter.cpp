/////////////////////////////////////////
// Converter.cpp
// simple calculator for min MUM length
/////////////////////////////////////////

// See the LICENSE file included with this software for license information.

//Converter.cpp
#include "Converter.h"  
using namespace std;   
int Converter( string infix, string &output, int size)  
{     
	Stack<char> convert;  
	//instantiation of char Stack convert     
	char b = 'b';      
	int c = 0;       
	infix = infix + ")";     
	convert.push('(');     
	int  i = 0;     
	int  len = infix.length();    
	while(!(convert.empty()) && i < len)    
	{      
		c = 0;  
		//reset '(' flag       
		//char to compare cases      
		switch(infix[i])      
		{        
			//if '(', place it onto the stack      
		case '(':        
			convert.push('(');        
			break;      
		case ')':        
			while(c != 1)//testing for '('        
			{ 	 
				b = convert.pop(); 	 
				if( b != '(') 
					//if not '(', output 	   
					output = output + b; 	 
				else 	   
					c = 1; 
				//else note in c that element was '('        
			}        
			break;      

		case '+':        
			while(convert.peek() != '(')        
			{ 	 
				b = convert.pop(); 	 
				output = output + b;        
			}        
			convert.push('+');        
			break;         
			//if '-', first pop all other operators from stack      
			//then push it onto the stack      
		case '-':        
			if((i == 0) && (infix[0] == '-'))        
			{ 	 
				//cout << '-'; 	 
				output = output + "-"; 	 
				break;        
			}        
			while(convert.peek() != '(')        
			{ 	 
				//cout << convert.pop(); 	 
				b = convert.pop(); 	 
				output = output + b;        
			}        
			convert.push('-');        
			break;      
			//if '*', first pop all other operators, except      
			// + or -, from stack      
			//then push it onto the stack      
		case '*':        
			while( convert.peek() != '(' )        
			{	    	 
				if((convert.peek() != '+')&&(convert.peek() != '-'))  
				{ 	   
					//cout << convert.pop(); 	   
					b = convert.pop(); 	   
					output = output + b; 	 
				} 	 
				else 	   
					break;        
			}        
			convert.push('*');        
			break;      
			//if '/', first pop all other operators, except      
			// + or -, from stack      
			//then push it onto the stack      
		case '/':        
			while(convert.peek() != '(')        
			{    	 
				if((convert.peek() != '+')&&(convert.peek() != '-')) 
				{ 	   
					//cout << convert.pop(); 	   
					b = convert.pop(); 	   
					output = output + b; 	 
				} 	 
				else 	   
					break;        
			}        
			convert.push('/');        
			break;      
		case '^':        
			while(convert.peek() != '(')        
			{  	 
				if(convert.peek() == '^') 	 
				{ 	   
					b = convert.pop(); 	   
					output = output + b; 	 
				} 	 
				else 	   
					break;         
			}        
			convert.push('^');        
			break;       
		case 'L':        
			while(convert.peek() != '(')        
			{  	 
				convert.pop(); 	 
				convert.pop();       
			}        
     
			convert.push('L');           
			break;      
			//ignore tabs, spaces from input	 
		case '\t':    
		case  ' ':      
			break;      

		default:      
			if(!(convert.empty()))      
			{  
				if(isdigit(infix[i]) || isalpha(infix[i]) ) 	 
				{ 	 
					output = output + infix[i]; 	  
					if(infix[i+1] != '.' && !(isdigit(infix[i+1]))) 	 
					{ 	     
						output = output + " "; 	   
					} 	 
				}      
			}       
			//to allow for floats      
			if(infix[i] == '.')    
			{ 	 	 
				output = output + ".";    
			}      
			break;   
   }     
   i++;   
 }   
 //cout << output << endl;    
 return 0;  
}    
float Calculator(string input, int size,float seqlen)  
{     
	Stack<float> calc; 
	//Stack calc of type float     
	char a;  
	//a to hold 1 char of input      
	float x, y, z;      
	int checkfornegative = 0;     
	int i = 0;      
	//while (cin >> a)//while there is still input to be read     
	while (i<=size)//while there is still input to be read     
	{       
		a = input[i];       
		//this is necessary if the first       
		//number in the expression is negative      
		//otherwise too many pops would be       
		//performed       
		while(checkfornegative == 0)       
		{         
			if(a == '-')         
			{ 	  
				--i; 	  
				float n; 	  
				n = input[i++]; 	 
				//cin >> n; 	  
				calc.push(n);       
			}        
			checkfornegative = 1;      
			//i++;    
		}       
		//determine if the first char is a digit     
		//or if it is a '.' to allow for use of floats    
		//in expression      
		if((isdigit(a)) || (a == '.'))       
		{        
			//if so, putback the digit and store the     
			//whole number into float b        
			//then push it into the stack      
			float b;        
			string s;      
			int start = i;       
			//cin.putback(a);     
			while ( isdigit (input[i] ) )      
			{ 	 //cout << input[i] << endl; 	 
				i++; 	  
				if ( input[i] == '.' ) 	  
					i++;       
			}       
			s = input.substr(start,i-start); 
			//cin >> b;         
			//cout << s << endl;      
			b = atof(s.c_str());        
			calc.push(b);     
		}      
		else //if it is not a digit     
			//figure out which operator it is      
 
			switch(a)        
			{         
			case 'S':        
			case 's': 	  
				calc.push(seqlen); 	 
				break;        
			case 'N':        
			case 'n': 	  
 
				break;         
			case 'M':         
			case 'm': 	  
				break;         
			case '+': 	 
				x = calc.pop(); 
				//right hand of expression 	 
				y = calc.pop(); 	 
				z = y + x;  
				calc.push(z); 
				//push value back onto stack  
				break;        
			case '-': 	 
				x = calc.pop();  
				y = calc.pop(); 	 
				z = y - x; 	 
				calc.push(z);        
				break;         
			case '*': 	  
				x = calc.pop();  
				y = calc.pop(); 	 
				z = y * x;  
				calc.push(z); 	 
				break;         
			case '/': 	 
				x = calc.pop(); 	 
				y = calc.pop(); 	 
				//testing for division by 0  
				if(x == 0) 	 
				{ 	    
					//cout << "Cant divide by 0!" << endl; 	   
					exit(1); 	 
				}          
				else 	  
				{ 	    
					z = y / x; 	    
					calc.push(z); 	 
				}          
				break;         
			case '^': 	 
				x = calc.pop();  
				y = calc.pop(); 	  
				z = pow(y,x); 	 
				calc.push(z); 	 
				break;         
			case 'L': 	 

				x = calc.pop(); 	 
				z = log(x)/log(2.0);  
				calc.push(z); 	  
				break;        
				//ignore tabs and spaces      
			case '\t':        
			case ' ':      
			case 'o':      
			case 'g':  
				break;    
			default: 	 
				break;       
  } 
  //end of switch statement       
  i++;     
}    
float out = calc.pop();   
out = ceil(out);     
return out; }  

