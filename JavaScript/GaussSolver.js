//GitHub: HenriqueIni
//www.blogcyberini.com

/*
Algoritmo para resolução de sistemas lineares via eliminação de Gauss
Complexidade no tempo: O(n^3)
A é a matriz dos coeficientes do sistema
b é a matriz dos coeficientes dos termos independentes
Forma do sistema (matricial): Ax = b
*/
function gaussSolver(A, b){
    var i, j, k, l, m;
    //ETAPA DE ESCALONAMENTO
    for(k = 0; k < A.length - 1; k++){
        //procura o maior k-ésimo coeficiente em módulo
        var max = Math.abs(A[k][k]);
        var maxIndex = k;
        for(i = k + 1; i < A.length; i++){
            if(max < Math.abs(A[i][k])){
                max = Math.abs(A[i][k]);
                maxIndex = i;
            }
        }
        if(maxIndex != k){
            /*
             troca a equação k pela equação com o
             maior k-ésimo coeficiente em módulo
             */
            for(j = 0; j < A.length; j++){
                var temp = A[k][j];
                A[k][j] = A[maxIndex][j];
                A[maxIndex][j] = temp;
            }
            var temp = b[k];
            b[k] = b[maxIndex];
            b[maxIndex] = temp;
        }
        //Se A[k][k] é zero, então a matriz dos coeficiente é singular
        //det A = 0
        if(A[k][k] == 0){
            return null;
        }else{
            //realiza o escalonamento
            for(m = k + 1; m < A.length; m++){
                var F = -A[m][k] / A[k][k];
                A[m][k] = 0; //evita uma iteração
                b[m] = b[m] + F * b[k];
                for(l = k + 1; l < A.length; l++){
                    A[m][l] = A[m][l] + F * A[k][l];
                }
            }
        }
    }
    //ETAPA DE RESOLUÇÃO DO SISTEMA
    var X = [];
    for(i = A.length - 1; i >= 0; i--){
        X[i] = b[i];
        for(j = i + 1; j < A.length; j++){
            X[i] = X[i] - X[j] * A[i][j];
        }
        X[i] = X[i] / A[i][i];
    }
    return X;
}
