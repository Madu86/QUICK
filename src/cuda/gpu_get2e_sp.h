
__device__ __forceinline__ void addint_oshell(QUICKULL* oULL, QUICKULL* obULL,QUICKDouble Y, int III, int JJJ, int KKK, int LLL,QUICKDouble hybrid_coeff,  QUICKDouble* dense, QUICKDouble* denseb,int nbasis);

__device__ __forceinline__ void addint(QUICKULL* oULL, QUICKDouble Y, int III, int JJJ, int KKK, int LLL,QUICKDouble hybrid_coeff,  QUICKDouble* dense, int nbasis);

__device__ __forceinline__ QUICKDouble quick_dsqr(QUICKDouble a);

__device__ __forceinline__ void FmT(int MaxM, QUICKDouble X, QUICKDouble* YVerticalTemp);

__device__ __forceinline__ int lefthrr(QUICKDouble RAx, QUICKDouble RAy, QUICKDouble RAz,
                                       QUICKDouble RBx, QUICKDouble RBy, QUICKDouble RBz,
                                       int KLMNAx, int KLMNAy, int KLMNAz,
                                       int KLMNBx, int KLMNBy, int KLMNBz,
                                       int IJTYPE,QUICKDouble* coefAngularL, int* angularL);

__device__ __forceinline__ int lefthrr23(QUICKDouble RAx, QUICKDouble RAy, QUICKDouble RAz,
                                    QUICKDouble RBx, QUICKDouble RBy, QUICKDouble RBz,
                                    int KLMNAx, int KLMNAy, int KLMNAz,
                                    int KLMNBx, int KLMNBy, int KLMNBz,
                                    int IJTYPE,QUICKDouble* coefAngularL, int* angularL);

__device__ __forceinline__ int lefthrr23_new(QUICKDouble RAx, QUICKDouble RAy, QUICKDouble RAz,
                                         QUICKDouble RBx, QUICKDouble RBy, QUICKDouble RBz,
                                         int KLMNAx, int KLMNAy, int KLMNAz,
                                         int KLMNBx, int KLMNBy, int KLMNBz,
                                         int IJTYPE,QUICKDouble* coefAngularL, int* angularL);
