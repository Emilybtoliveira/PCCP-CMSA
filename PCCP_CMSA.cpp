#include <bits/stdc++.h>
#include <ilcplex/ilocplex.h>
#include <chrono>
using namespace std;

int n, m, t, min_alpha, max_alpha; // numero de cadeias, tamanho de cada cadeia, tamanho do alfabeto
//========================== estruturas imutaveis ===========================
vector<char> dataset_alphabet;
vector<string> strings_dataset;
vector<vector<int>> integer_dataset;
map<char, int> alphabetMap;

//========================== estruturas mutaveis a cada loop ===========================
vector<vector<int>> columns_sets;
vector<int> columns_sets_ages;
vector<vector<int>> sets_using_columns;
vector<vector<int>> columns_sets_ham_distances;
vector<int> columns_ILP_selection;
vector<vector<int>> sets_closest_strings;

typedef IloArray<IloNumVarArray> NumVar2D;
typedef IloArray<NumVar2D> NumVar3D;
#define INF 2147483647;

//========================DEBUG FUNCTIONS: DELETE LATER===================//

void printSets(vector<int> solution)
{
    cout << "colunas por conjunto:\n";
    for (vector<int> h: columns_sets)
    {
        for(int i: h){
            cout << i << " ";
        }
        cout << endl;
    }  

    cout << "conjuntos por coluna\n";
    for (vector<int> g: sets_using_columns)
    {
        for(int i: g){
            cout << i << " ";
        }
        cout << endl;
    } 
/* 
    for(int i = 0; i < (int)solution.size(); i++){        
        cout << "indice " << i << ": " << solution[i] << endl;
        if(solution[i] == 1){
            for (int j: columns_sets[i]){ cout << j << " ";}
            cout << endl;
        }
    }
    cout << endl; */
}

void mergeSetsHammingDistances(vector<int> solution, vector<vector<int>> distances)
{
    int d, maior_distancia = 0;
    vector<int> selected_sets;

    cout << "conjuntos selecionados:" << endl;

    for(int i = 0; i < (int)solution.size(); i++)
    {               
        if(solution[i] == 1)
        {
            cout << i << " ";
            selected_sets.push_back(i);
        }
    }
    cout << endl;

    cout << "distancias:" << endl;
    for(int i = 0; i < n; i++){
        d = 0;
        cout << i << ": ";
        for(int k: selected_sets)
        {
            d += columns_sets_ham_distances[k][i];
            cout << columns_sets_ham_distances[k][i] << " ";
        }
        
        cout << "= " << d << endl;

        if(d > maior_distancia)
        {
            maior_distancia = d;
        }
    }

    cout << maior_distancia << endl;
}

void getHammingDistance(vector<int> closest_string)
{
    long unsigned int i;
    int d = 0;
    int maior_distancia = 0;

    for (vector<int> instance: integer_dataset)
    {
        d = 0;

        for (i = 0; i < closest_string.size(); i++)
        {   
           // cout << "instancia: " << instance[i]<< " solução: " << closest_string[i] << endl;

            if (instance[i] != closest_string[i])
            { 
                d += 1; 
            } 
        }
        
      //  cout << "distancia: " << d << endl << endl << endl;

        if(d > maior_distancia){
            maior_distancia = d;
        }
    }

    cout << maior_distancia << endl;
}

void generateMergedClosestString(vector<int> solution){
    vector<int> selected_sets;
    vector<int> closest_string;
    
    for(int i = 0; i < (int)solution.size(); i++)
    {               
        if(solution[i] == 1)
        {            
            selected_sets.push_back(i);
        }
    }

    for(int i: selected_sets)
    {
        cout << "(";
        for(int j: columns_sets[i])
        {
            cout << j << ",";
        }
        cout << ") = ";

        cout << "[";
        for(int j: sets_closest_strings[i])
        {
            cout << j << ",";
        }
        cout << "]"<< endl;
    }
    

    for (int i = 0; i < m; i++)  
    {
        closest_string.push_back(0);   
    }

    for(int i: selected_sets)
    {
        for(int j = 0; j < (int)columns_sets[i].size(); j++)
        {
                //cout << columns_sets[i][j] << " " << sets_closest_strings[i][j] << endl;

                closest_string[columns_sets[i][j]] = sets_closest_strings[i][j];
        }
    } 

    cout << "closest string: ";
    for(int i: closest_string){
        cout << i << " ";
    }
    cout << endl;

    getHammingDistance(closest_string);

}

//======================== CMSA CORE FUNCTIONS ===================//

void initializeDS()
{
    columns_ILP_selection.clear();

    for (long unsigned int i = 0; i < columns_sets.size(); i++)  
    {
        columns_ILP_selection.push_back(0);   
    }
}

int setsSelectionSolver(int n_sets, vector<vector<int>> C, vector<vector<int>> e)
{
    IloEnv env;
    IloModel Model(env, "Problema da Seleção de Colunas");
    int cost = 0;
    
    try {      
        //Variável de decisão
		IloIntVarArray x(env, n_sets, 0, 1);

        //Função objetivo.
		IloNumVar z(env, 0, IloInfinity, ILOINT);
        IloExpr exp3(env);
        exp3 = z;
		Model.add(IloMinimize(env, exp3));

        //Restrições
        for (int i = 0; i < m; i++){
            IloExpr exp1(env);
            
            //cout << C[i].size() << endl;   
            for(long unsigned int s = 0; s < C[i].size(); s++){               
                exp1 += x[C[i][s]];
                //cout << exp1;
            }

            Model.add(exp1 == 1);
        }

        for (int j = 0; j < n; j++){
            IloExpr exp2(env);

            for (int s = 0; s < n_sets; s++){
                exp2 += (e[s][j] * x[s]);
            }

            Model.add(exp2 <= z);
        }

        //Solving
        IloCplex cplex(Model);
		cplex.setOut(env.getNullStream());
       // cout << cplex.getModel() << endl;
		
		if (!cplex.solve()) {
			env.error() << "Falhou ao otimizar o problema" << endl;
		}

        //cout << cplex.getValue(z) << endl ;
        cost = cplex.getObjValue();
		//cout << "custo otimo: " << cost << endl; 
		
		//Obtendo a solução
		IloNumArray sol(env, n_sets);
		cplex.getValues(sol, x);

        initializeDS();

		for (int i = 0; i < n_sets; i++){
            columns_ILP_selection[i] = sol[i];
            /* if (((int)sol[i]) >= 1){
                columns_ILP_selection[i] = 1;
            }
            else{
                columns_ILP_selection[i] = 0;
            } */
			cout << columns_ILP_selection[i] << " ";			
        }        
        cout << endl;

    } catch (const IloException& e) {
		cerr << "Exception caught: " << e << endl;
	} catch (...) {
		cerr << "Unknown exception caught!" << endl;
	}

    env.end();
    
    //==================DEBUG FUNCTIONS=============================
   // printSets(columns_ILP_selection);
    generateMergedClosestString(columns_ILP_selection);
    
    return cost;
}

vector<int> PCCPSolver(int n, int m, int min_alpha, int max_alpha, vector<vector<int>> S)
{
    int delta = max_alpha - min_alpha;
	vector<int> closest_string;
    
    // Criando o ambiente
	IloEnv env;
	IloModel Model(env, "Problema da Cadeia");

	try
	{
// region Variável de decisão

		NumVar2D z(env, n);
		IloNumVarArray t(env, m, min_alpha, max_alpha, ILOINT);
		IloNumVar d(env, 0, IloInfinity, ILOINT);

		for (int i = 0; i < n; i++) {
			z[i] = IloNumVarArray(env, m, 0, 1, ILOINT);
		}


// region Função Objetivo
		
		IloExpr exp0(env);
		exp0 = d;		
		Model.add(IloMinimize(env, exp0));


// region Restrições

		for (int i = 0; i < n; i++)
		{
			IloExpr exp1(env);
			for (int j = 0; j < m; j++)
			{
				exp1 += z[i][j];
			}

			Model.add(exp1 <= d);
		}
		
		for (int i = 0; i < n; i++)
		{
			IloExpr exp2(env);
			for (int j = 0; j < m; j++)
			{
				exp2 = t[j] - S[i][j];
				Model.add(exp2 <= delta*z[i][j]);
			}		
		}

		for (int i = 0; i < n; i++)
		{
			IloExpr exp3(env);
			for (int j = 0; j < m; j++)
			{
				exp3 = S[i][j] - t[j];
				Model.add(exp3 <= delta * z[i][j]);
			}
		}

		IloCplex cplex(Model);
		cplex.setOut(env.getNullStream());

		
		if (!cplex.solve()) {
			env.error() << "Falhou ao otimizar o problema" << endl;
		}

		/* double obj = cplex.getObjValue();
		cout << "custo otimo: " << obj << endl;  */
		
		 //Obtendo a solu��o
		IloNumArray sol(env, m);
		cplex.getValues(sol, t);
        

		// Imprimindo a solu��o
		for (int i = 0; i < m; i++)
        {
            if (sol[i] < 1){
                sol[i] = 0;
            }			
            //cout << sol[i] << " ";            
            closest_string.push_back((int)sol[i]);
        }
        //cout << endl; 
	}
	catch (const IloException& e)
	{
		cerr << "Exception caught: " << e << endl;
	}
	catch (...)
	{
		cerr << "Unknown exception caught!" << endl;
	} 
    

	env.end();
    return closest_string;
}

void printAges(){
    cout << "ages:\n";
    
    for (int h: columns_sets_ages)
    {
        cout << h << " ";
    }  
    cout << endl;
    cout << "fim ages\n";

}

void columnsSelectorv2(int loops, int max_size)
{
    srand(time(0));
    int set, i, j, l;
    int root_index = columns_sets.size();

    float columns_partition_size = (rand() % max_size) + 1;
    int n_sets_per_loop = ceil((float) m / columns_partition_size);
    int total_n_sets = n_sets_per_loop * loops;
  
    vector<int> shuffled_sets;
    vector<int> empty_vec;   
    
    
    //================= pré processamento ==========================
    // cout << "conjuntos por loop: " << n_sets_per_loop << endl;
     
    //inicializando o vetor de colunas pra poder acessar via indice
    for (i = 0; i < total_n_sets; i++)
    {
        columns_sets.push_back(empty_vec);
        columns_sets_ages.push_back(0);   
    }
    
    //================= loop da criação dos conjuntos ==========================
    for (l = 0; l < loops; l++)
    {        
        int starting_index = (l* n_sets_per_loop) + root_index;

        //popula o shuffled_sets com os indices dos conjuntos em quantidade igual a columns_partition_size        
        for(i = starting_index; i < (starting_index + n_sets_per_loop); i++)
        {
            //seta as idades dos novos conjuntos para 0
            for(j = 0; j < columns_partition_size; j++)
            {
                shuffled_sets.push_back(i);
            }
        }

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        shuffle (shuffled_sets.begin(), shuffled_sets.end(), std::default_random_engine(seed));
        
        /* cout << "conjuntos shuffled\n";
        for (int i: shuffled_sets)
        {
            cout << i << " ";
        }
        cout << endl; */

        
        /* aqui cada coluna é colocada em um vetor de columns_sets na ordem em que eles aparecem em shuffled_sets, garantindo a aleatoriedade 
            das colunas nos conjuntos, que cada conjunto nao ultrapassa seu tamanho e que as colunas serao inseridas em ordem crescente */
        j = 0;            
        for (i = 0; i < m; i++)
        {
            //column = rand() % n_sets_per_loop;       
            set = shuffled_sets[j];
            columns_sets[set].push_back(i);

            //sets_using_columns[i].push_back(set);  
            j++;  
        }

        shuffled_sets.clear(); 
        //cout << endl;       
    }
   
   /*  cout << "colunas por conjunto:\n";
    for (vector<int> h: columns_sets)
    {
        for(int i: h){
            cout << i << " ";
        }
        cout << endl;
    }  
    
     */
   // printAges();
}

void columnsSelector(float columns_partition_size, int loops, int root_index)
{
    int set, i, j, l;
    int n_sets_per_loop = ceil((float) m / columns_partition_size);
    int total_n_sets = n_sets_per_loop * loops;
  
    vector<int> shuffled_sets;
    vector<int> empty_vec;    
    srand(time(0));
    
    
    //================= pré processamento ==========================
    // cout << "conjuntos por loop: " << n_sets_per_loop << endl;
     
    //inicializando o vetor de colunas pra poder acessar via indice
    for (i = 0; i < total_n_sets; i++)
    {
        columns_sets.push_back(empty_vec);
        columns_sets_ages.push_back(0);   
    }
    
    //================= loop da criação dos conjuntos ==========================
    for (l = 0; l < loops; l++)
    {        
        int starting_index = (l* n_sets_per_loop) + root_index;

        //popula o shuffled_sets com os indices dos conjuntos em quantidade igual a columns_partition_size        
        for(i = starting_index; i < (starting_index + n_sets_per_loop); i++)
        {
            //seta as idades dos novos conjuntos para 0
            for(j = 0; j < columns_partition_size; j++)
            {
                shuffled_sets.push_back(i);
            }
        }

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        shuffle (shuffled_sets.begin(), shuffled_sets.end(), std::default_random_engine(seed));
        
        /* cout << "conjuntos shuffled\n";
        for (int i: shuffled_sets)
        {
            cout << i << " ";
        }
        cout << endl; */

        
        /* aqui cada coluna é colocada em um vetor de columns_sets na ordem em que eles aparecem em shuffled_sets, garantindo a aleatoriedade 
            das colunas nos conjuntos, que cada conjunto nao ultrapassa seu tamanho e que as colunas serao inseridas em ordem crescente */
        j = 0;            
        for (i = 0; i < m; i++)
        {
            //column = rand() % n_sets_per_loop;       
            set = shuffled_sets[j];
            columns_sets[set].push_back(i);

            //sets_using_columns[i].push_back(set);  
            j++;  
        }

        shuffled_sets.clear(); 
        //cout << endl;       
    }
   
   /*  cout << "colunas por conjunto:\n";
    for (vector<int> h: columns_sets)
    {
        for(int i: h){
            cout << i << " ";
        }
        cout << endl;
    }  
    
    printAges(); */
}

void createSetsUsingColumns()
{
    //limpa o conjunto
    sets_using_columns.clear();

    //inicializando o vetor de conjuntos por coluna pra poder acessar via indice
    vector<int> empty_vec;    
    for (int i = 0; i < m; i++)  
    {
        sets_using_columns.push_back(empty_vec);   
    }

    //atribui conjuntos às colunas
    for (long unsigned int i = 0; i < columns_sets.size(); i++)
    {
        for(int j: columns_sets[i]){
            sets_using_columns[j].push_back(i);  
        }
    } 

   /*  cout << "conjuntos por coluna\n";
    for (vector<int> g: sets_using_columns)
    {
        for(int i: g){
            cout << i << " ";
        }
        cout << endl;
    }   */
}

vector<vector<int>> constructReducedInstance(vector<int> columns_set)
{
    vector<vector<int>> reduced_instance;
    vector<int> instance;

    for (int i = 0; i < n; i++)
    {
        for (int j: columns_set)
        {
            instance.push_back(integer_dataset[i][j]);
        }

        reduced_instance.push_back(instance);
        instance.clear();
    }

/*     for (vector<int> instance: reduced_instance)
    {
        for (int integer: instance)
        {
            cout << integer << " ";
        }
        cout << endl;
    } */

    return reduced_instance;
}

void computeHammingDistance(vector<vector<int>> reduced_instance, vector<int> closest_string)
{
    long unsigned int i;
    int d = 0;
    vector<int> ham_distance;

    for (vector<int> instance: reduced_instance)
    {
        d = 0;

        for (i = 0; i < closest_string.size(); i++)
        {   
            if (instance[i] != closest_string[i])
            { 
                d += 1; 
            } 
        }
        
        //cout << d << " ";

        ham_distance.push_back(d);
    }
   // cout << endl;

    columns_sets_ham_distances.push_back(ham_distance);
}

void solveSmallInstances()
{
    vector<int> empty_vec;
    vector<vector<int>> reduced_instance;
    vector<int> closest_string;
    
    for(vector<int> set: columns_sets)
    {
        reduced_instance = constructReducedInstance(set);
        closest_string = PCCPSolver(n, set.size(), min_alpha, max_alpha, reduced_instance);
        sets_closest_strings.push_back(closest_string);
        computeHammingDistance(reduced_instance, closest_string);

        reduced_instance.clear();
        closest_string.clear();
        //return;
    }


   /*  cout << "distancias:\n";
    for(vector<int> set_distance: columns_sets_ham_distances)
    {
        for(int distance: set_distance)
        {
            cout << distance << " "; 
        }
        cout << endl;
    } */
}

void adapt(int max_age)
{
    auto iterator_1 = columns_sets.begin();
    auto iterator_2 = columns_sets_ages.begin();

    vector<int> index_remove;

    for (long unsigned int i = 0; i < columns_sets_ages.size(); i++)
    {
        if(columns_ILP_selection[i] == 0)
        {
            columns_sets_ages[i] += 1;

            if (columns_sets_ages[i] > max_age){
                //cout << "deletando o " << i << "com idade " << columns_sets_ages[i] << endl;
                index_remove.push_back(i);
            }
        } 
        else 
        {
            columns_sets_ages[i] = 0;
        }
    }

    for (int i = 0; i < (int)index_remove.size(); i++)
    {
        //int index = ()
        columns_sets.erase(iterator_1 +  index_remove[i] - i);
        columns_sets_ages.erase(iterator_2 + index_remove[i] - i);
    }

    //printAges();
}

void generateAlphabetMapping()
{
    int i;
    for (i = 0; i < t; i++)
    {
        char alpha_char = dataset_alphabet[i];
        alphabetMap[alpha_char] = i;
    }

    min_alpha = 0;
    max_alpha = i-1;

    /* for (const auto &item : alphabetMap)
    {
        cout << "[" << item.first << ", " << item.second << "]\n";
    } */
}

void instanceTransformFunc()
{
    vector<int> int_string;

    for (string i: strings_dataset)
    {
        int_string.clear();

        for (char c: i) 
        {
            int_string.push_back(alphabetMap[c]);
        }
        integer_dataset.push_back(int_string);
    }

    /* for(vector<int> j: integer_dataset){
        for(int inteiro: j){
            cout << inteiro << " ";
        }
        cout << "\n\n";
    } */
}

void readInstance()
{
    cin >> n >> m >> t; // numero de cadeias, tamanho das cadeias e tamanho do alfabeto

    char cur_char;

    for (int i = 0; i < t; i++)
    {
        cin >> cur_char;
        dataset_alphabet.push_back(cur_char);
    }

    for (int i = 0; i < n; i++)
    {
        string cur_string;
        cin >> cur_string;
        strings_dataset.push_back(cur_string);
    }
}

void printInstance()
{
    printf("%d %d %d\n", n, m, t);

    for (char i : dataset_alphabet)
    {
        cout << i << '\n';
    }

    for (string i : strings_dataset)
    {
        cout << i << '\n';
    }
}

void printSetsSolution(vector<int> solution){

    for(int i = 0; i < (int)solution.size(); i++){        
        cout << "indice " << i << ": " << solution[i] << endl;
        if(solution[i] == 1){
            for (int j: columns_sets[i]){ cout << j << " ";}
            cout << endl;
        }
    }
    cout << endl;

    cout << "conjuntos por coluna\n";
    for (vector<int> g: sets_using_columns)
    {
        for(int i: g){
            cout << i << " ";
        }
        cout << endl;
    }  

    /*  cout << "colunas por conjunto:\n";
    for (vector<int> h: columns_sets)
    {
        for(int i: h){
            cout << i << " ";
        }
        cout << endl;
    }   */
       
}

void mainLoop()
{
   /*  //se a instância for pequena, a partição não vai dar certo
    if (m <= 10)
    { 
        PCCPSolver(n, m, min_alpha, max_alpha, integer_dataset);
        return;
    } */
    
    int opt, bsf = INF;
    int loops = 0;
    vector<int> bsf_selection;

    auto loop_start = chrono::high_resolution_clock::now();
    auto loop_end = loop_start;
    auto loop_cur = chrono::duration_cast<chrono::milliseconds>(loop_end - loop_start);

    while(loop_cur.count() < 60000)
    //while(loops < 10)
    {
        //================= CMSA Procedures ==========================
        //columnsSelector(2, 1, columns_sets.size());
        columnsSelectorv2(1, 6);
        createSetsUsingColumns();            
        solveSmallInstances();
        
        opt = setsSelectionSolver(columns_sets.size(), sets_using_columns, columns_sets_ham_distances); 
        cout << bsf << ", " << opt << endl;

        //printSetsSolution(columns_ILP_selection);

        if(opt < bsf){
            bsf = opt;
            bsf_selection = columns_ILP_selection;
        }

        adapt(4);

        //limpa estruturas de dados que vao ser reutilizadas nos loops
        columns_sets_ham_distances.clear();
        columns_ILP_selection.clear();
        sets_closest_strings.clear();
        
        loops += 1;

        loop_end = chrono::high_resolution_clock::now();
        loop_cur = chrono::duration_cast<chrono::milliseconds>(loop_end - loop_start);
    }
   // printAges();
   //cout << "loops: " << loops << endl; 
   //printSetsSolution(bsf_selection);
}

int main()
{
    readInstance();
    //printInstance();
    generateAlphabetMapping();
    instanceTransformFunc();
    
    /* vector<int> t = PCCPSolver(n, m, min_alpha, max_alpha, integer_dataset);
    getHammingDistance(t);
    return 0; */

    auto start = chrono::high_resolution_clock::now(); 
    mainLoop();
    auto end = chrono::high_resolution_clock::now();
    auto Elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);


    cout << "Tempo decorrido(ms): " << Elapsed.count() << endl;
}
