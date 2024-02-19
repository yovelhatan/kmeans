#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**************** Function declarations ****************/

typedef struct{
    double ** points;
    int size;
} Cluster;

void malloc_failed();
int error_handler(int K, int N, int d, int iter);
double ** read_points_from_stdin(int N, int d);
double ** init_centroids(int d, int K, double ** points_array);
int argmin_from_distances(int K, int d, double * point, double ** centroids);
Cluster ** add_points_to_clusters(double ** points_array, double ** centroids, int N, int d, int K);
double get_distance(double point[], double centroid[], int d);
double ** update_centroids(int K, int d, Cluster ** clusters);
int is_delta_centroids_small(double ** new_centroids, double ** old_centroids, int K, int d);
void print_centroids(double ** centroids, int K, int d);

/**************** Main *********************************/

int main(int argc, char * argv[]) {
    const int MAX_ITER = 200;
    int K, N, d, iter;
    double ** points_array;
    double ** centroids;
    double ** new_centroids;
    Cluster ** clusters;

    if (argc != 4 && argc != 5) {
        printf("An Error Has Occurred");
        exit(1);
    }

    K = atoi(argv[1]);
    N = atoi(argv[2]);
    d = atoi(argv[3]);
    iter = (argc == 5) ? atoi(argv[4]) : MAX_ITER;
    
    if (error_handler(K, N, d, iter) == 1){
        exit(1);
    }

    points_array = read_points_from_stdin(N, d);
    centroids = init_centroids(d, K, points_array);
    clusters = NULL;
    new_centroids = NULL;

    while (iter) {
        clusters = add_points_to_clusters(points_array, centroids, N, d, K);
        new_centroids = update_centroids(K, d, clusters);

        if (is_delta_centroids_small(new_centroids, centroids, K, d) == 0) {
            break;
        }
        
        centroids = new_centroids;
        iter -= 1;
    }
    
    print_centroids(centroids, K, d);
    
    free(points_array);
    free(clusters);
    free(new_centroids);
    free(centroids);

    return 0;
}

/**************** Implementaions ***********************/

/*
Memory allocation failed
*/
void malloc_failed() {
    fprintf(stderr, "Memory allocation failed\n");
    exit(1); 
}

/*
Insert each point coordinates array to points array
*/
double **read_points_from_stdin(int N, int d) {
    int i, j;

    /* Allocate memory for points array */
    double ** points_array = (double **)malloc(N * sizeof(double *));
    double * coordinates;

    /* Check if the allocation was successful */
    if (points_array == NULL) {
        malloc_failed();
    }

    coordinates = (double *)malloc(N * d * sizeof(double));

    /* Check if the allocation was successful */
    if (coordinates == NULL) {
        free(points_array);
        malloc_failed(); 
    }

    for (i = 0; i < N; i++) {
        points_array[i] = coordinates + i * d;
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j < d; j++) {
            if (j == d - 1){
                if (scanf("%lf", &points_array[i][j]) != 1) {
                    free(coordinates); 
                    free(points_array);
                    malloc_failed();
                }
            }
            else if (scanf("%lf,", &points_array[i][j]) != 1) {
                free(coordinates); 
                free(points_array);
                malloc_failed(); 
            }
        }
    }

    return points_array;
}

/*
Check the parameters are valid 
*/
int error_handler(int K, int N, int d, int iter){
    int failed = 0;
  
    if ((int)K != K  || K <= 1 || K >= N){
        fprintf(stderr, "Invalid number of clusters!\n");
        failed = 1;
    }

    if ((int)N != N  || N <= 1 ){
        fprintf(stderr,"Invalid number of points!\n");
        failed = 1;
    }

    if ((int)d != d){
        fprintf(stderr, "Invalid dimension of point!\n");
        failed = 1;
    }
    
    if ((int)iter != iter || iter <= 1 || iter >= 1000){
        fprintf(stderr, "Invalid maximum iteration!\n");
        failed = 1;
    }

    return failed;
}

/*
Get the distance between a point and a centroid
*/
double get_distance(double point[], double centroid[], int d) {
    double distance = 0;
    int i;

    for(i = 0; i < d; i++) {
        distance += pow(point[i] - centroid[i], 2);
    }

    return sqrt(distance);
}

/*
Initialize the first centroids
*/
double ** init_centroids(int d, int K, double ** points_array) {
    double ** centroids = (double **)malloc(K * sizeof(double *));
    int i, j;
    double * coordinates;

    /* Check if the allocation was successful */
    if (centroids == NULL) {
        malloc_failed();
    }

    coordinates = (double *)malloc(K * d * sizeof(double));

    /* Check if the allocation was successful */
    if (coordinates == NULL) {
        free(centroids);
        malloc_failed(); 
    }
    
    for (i = 0; i < K; i++) {
        for (j = 0; j < d; j++) {
            centroids[i + j * d] = points_array[i + j * d];
        }
    }

    return centroids;
}

/*
Argmin function - return the index of the minimal distance from point
*/
int argmin_from_distances(int K, int d, double * point, double ** centroids) {
    int i, cluster_index = 0;
    double * distances = calloc(K, sizeof(double));

    for (i = 0; i < K; i++) {
        distances[i] = get_distance(point, centroids[i], d);
        if (distances[i] < distances[cluster_index]) {
            cluster_index = i;
        }
    }
    
    free(distances);
    return cluster_index;
}

/*
For all points check which cluster is right and add it to there
*/
Cluster ** add_points_to_clusters(double ** points_array, double ** centroids, int N, int d, int K) {
    Cluster ** clusters = (Cluster **)calloc(K, sizeof(Cluster *));
    int i, j;

    if (clusters == NULL){
        malloc_failed();
    }
    
    /* Allocate memory for each cluster */
    for (i = 0; i < K; i++) {
        clusters[i] = (Cluster *)calloc(1, sizeof(Cluster));

        if (clusters[i] == NULL){
            /* Free memory allocated so far */
            for (j = 0; j < i; j++) {
                free(clusters[j]->points);
                free(clusters[j]);
            }
            free(clusters);
            malloc_failed();
        }

        clusters[i]->size = 0;
    }

    /* Add each point to the right cluster */
    for (i = 0; i < N; i++) {
        int cluster_index = argmin_from_distances(K, d, points_array[i], centroids);
        int new_size = clusters[cluster_index]->size + 1;
        
        double ** temp_cluster_points = (double **)malloc(new_size * sizeof(double *));
        if (temp_cluster_points == NULL){
            malloc_failed();
        }

        for (j = 0; j < new_size - 1; j++) {
            temp_cluster_points[j] = clusters[cluster_index]->points[j];
        }

        temp_cluster_points[new_size - 1] = points_array[i];
        free(clusters[cluster_index]->points);
        clusters[cluster_index]->points = temp_cluster_points;
        clusters[cluster_index]->size = new_size;
    }
    
    return clusters;
}

/*
Update the centroids
*/
double ** update_centroids(int K, int d, Cluster ** clusters){
    /* Allocate memory for new_centroids array */
    double ** new_centroids = (double **)malloc(K * sizeof(double *));
    int i, j, k;

    /* Check if the allocation was successful */
    if (new_centroids == NULL){
        malloc_failed();
    }

    for (i = 0; i < K; i++){
        new_centroids[i] = (double *)malloc(d * sizeof(double));
        if (new_centroids[i] == NULL){
            /* Free memory allocated so far */
            for (j = 0; j < i; j++) {
                free(new_centroids[j]);
            }
            malloc_failed();
        }
    }

    /* Calculates new_centroid for each cluster */
    for (i = 0; i < K; i++){
        int number_of_points = clusters[i]->size;
        for (j = 0; j < d; j++){
            double sum = 0;
            for (k = 0; k < number_of_points; k++){
                sum += clusters[i]->points[k][j];
            }
            new_centroids[i][j] = sum / number_of_points;
        }
    }

    return new_centroids;
}

/*
Check new and old centroids
*/
int is_delta_centroids_small(double ** new_centroids, double ** old_centroids, int K, int d) {
    double EPSILON = 0.001;
    int i;

    for (i = 0; i < K; i++) {
        double distance = get_distance(new_centroids[i], old_centroids[i], d);
        if (distance > EPSILON) {
            return 1;
        }
    }

    return 0;
}

/*
Print end centroids with round to 4 digits after the dot
*/
void print_centroids(double ** centroids, int K, int d) {
    int i, j;
    
    for (i = 0; i < K; i++) {
        for (j = 0; j < d; j++) {
            printf("%.4f", centroids[i][j]);
            if (j < d - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
    printf("\n");
}
