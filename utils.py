import numpy as np
import matplotlib.pyplot as plt
import tensorly as tl
from tensorly import unfold
from tensorly import random
from tensorly.tenalg import multi_mode_dot, mode_dot
from tensorly.decomposition import tucker
from scipy.stats import norm
import json
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler

def run_tensor_regression_cv(region_id,     # region_id
                              tensor_path="expression_tensor.npz",
                              vcf_base_path="/projects/YangLabData/Ruilong/APE_project/plink_result/",
                              json_path="/projects/YangLabData/Ruilong/APE_project/valid_region_gene_map.json",
                              n_splits=5,
                              random_state=42):
    # === Load expression tensor ===
    loaded = np.load(tensor_path, allow_pickle=True)
    tensor_3d = loaded["tensor"]
    gene_list = loaded["gene_names"].tolist()
    unique_samples = loaded["sample_names"].tolist()
    cell_types = loaded["cell_types"].tolist()
    r2_scores = []
    y_true_list = []
    y_pred_list = []

    # === Slice tensor by region ===
    subtensor, subgenes = slice_tensor_by_region(region_id, tensor_3d, gene_list, json_path)
    subtensor = np.transpose(subtensor, (2, 1, 0))  # shape: (samples, cell_types, genes)

    # === Tucker ranks ===
    rank1 = cal_eigen_varexp(subtensor, mode0=1, modes=[0, 2], varexp=1)
    rank2 = cal_eigen_varexp(subtensor, mode0=2, modes=[0, 1], varexp=1)

    # === Read covariates ===
    vcf_path = f"{vcf_base_path}/{region_id}/{region_id}_final.raw"
    vcf_df = pd.read_csv(vcf_path, sep=r'\s+')
    vcf_df = vcf_df.fillna(vcf_df.mean(numeric_only=True))

    # === Construct design matrix ===
    X = vcf_df.set_index("FID").iloc[:, 5:]

    # === Standardize SNP covariates ===
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)

    # === Add intercept term ===
    X_scaled.insert(0, "intercept", 1)
    X_scaled.index.name = "Sample_ID"
    design_matrix = X_scaled

    n_samples, n_covariates = design_matrix.shape
    if n_samples <= n_covariates:
        raise ValueError(f"Design matrix invalid: {n_samples} samples < {n_covariates} covariates.")

    # === Cross-validation ===
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    r2_scores = []

    for fold, (train_idx, test_idx) in enumerate(kf.split(subtensor)):
        print(f"\n🔄 Fold {fold+1}/{n_splits}...")

        # Split data
        subtensor_train = subtensor[train_idx]
        subtensor_test = subtensor[test_idx]
        design_train = design_matrix.iloc[train_idx]
        design_test = design_matrix.iloc[test_idx]

        # Run tensor regression on training set
        core_shape = (design_train.shape[1], rank1, rank2)
        results = TensorReg_3d(
            tensor=subtensor_train,
            X_covar1=design_train,
            core_shape=core_shape
        )

        # Extract intermediate C_ts and test covariates
        C_ts = results["C_ts"]  # shape: (n_train, cell_types, genes)
        X_covar2 = np.eye(subtensor.shape[1])  # identity for cell types
        X_covar3 = np.eye(subtensor.shape[2])  # identity for genes

        # Predict test tensor using learned C_ts projection and test covariates
        U_test = multi_mode_dot(
            C_ts, [design_test.values, X_covar2, X_covar3], modes=[0, 1, 2]
        )

        y_true = subtensor_test.flatten()
        y_pred = U_test.flatten()

        r2 = r2_score(y_true, y_pred)
        r2_scores.append(r2)
        y_true_list.append(y_true)
        y_pred_list.append(y_pred)

    return r2_scores, y_true_list, y_pred_list




def run_tensor_regression(region_id,
                          tensor_path="expression_tensor.npz",
                          vcf_base_path="/projects/YangLabData/Ruilong/APE_project/plink_result/",
                          json_path="/projects/YangLabData/Ruilong/APE_project/valid_region_gene_map.json"):
    # === Load expression tensor ===
    loaded = np.load(tensor_path, allow_pickle=True)
    tensor_3d = loaded["tensor"]
    gene_list = loaded["gene_names"].tolist()
    unique_samples = loaded["sample_names"].tolist()
    cell_types = loaded["cell_types"].tolist()
    
    # === Slice tensor by region ===
    subtensor, subgenes = slice_tensor_by_region(region_id, tensor_3d, gene_list, json_path)
    
    # === Transpose to (samples, cell_types, genes) ===
    subtensor = np.transpose(subtensor, (2, 1, 0))

    # === Calculate Tucker ranks ===
    rank1 = cal_eigen_varexp(subtensor, mode0=1, modes=[0, 2], varexp=1)
    rank2 = cal_eigen_varexp(subtensor, mode0=2, modes=[0, 1], varexp=1)

    # === Read VCF-based covariates ===
    vcf_path = f"{vcf_base_path}/{region_id}/{region_id}_final.raw"
    vcf_df = pd.read_csv(vcf_path, sep=r'\s+')
    vcf_df = vcf_df.fillna(vcf_df.mean(numeric_only=True))
    
    # === Construct design matrix ===
    X = vcf_df.set_index("FID").iloc[:, 5:]

    # === Standardize SNP covariates ===
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)

    # === Add intercept term ===
    X_scaled.insert(0, "intercept", 1)
    X_scaled.index.name = "Sample_ID"
    design_matrix = X_scaled
    n_samples, n_covariates = design_matrix.shape
    if n_samples <= n_covariates:  
        raise ValueError(f"Design matrix invalid: {n_samples} samples < {n_covariates} covariates. Need more samples than covariates to perform regression.")
    
    # === Define core shape ===
    core_shape = (design_matrix.shape[1], rank1, rank2)

    # === Run tensor regression ===
    All_cov = TensorReg_3d(tensor=subtensor,
                           X_covar1=design_matrix,
                           core_shape=core_shape)
    
    return All_cov

def slice_tensor_by_region(region_id, tensor_3d, gene_list, genelist_json_path):
    """
    Slice a 3D gene expression tensor by region based on gene names.

    Parameters:
        region_id (str): e.g., "1_10583_1961168"
        tensor_3d (np.ndarray): shape = (num_genes, num_cell_types, num_samples)
        gene_list (list): gene names corresponding to axis 0 of tensor_3d
        genelist_json_path (str): path to JSON file mapping region_id → gene list

    Returns:
        sliced_tensor (np.ndarray): shape = (num_genes_in_region, num_cell_types, num_samples)
        selected_genes (list): the gene names found and used from the region
    """
    # Load the gene list JSON file
    with open(genelist_json_path, 'r') as f:
        genelist = json.load(f)

    if region_id not in genelist:
        raise ValueError(f"❌ region ID '{region_id}' not found in JSON file.")

    region_genes = genelist[region_id]['genes']

    # Map gene names to their indices in the tensor
    gene_to_index = {gene: idx for idx, gene in enumerate(gene_list)}

    # Find genes in the region that exist in the tensor
    selected_genes = [gene for gene in region_genes if gene in gene_to_index]
    selected_indices = [gene_to_index[gene] for gene in selected_genes]

    # Slice the tensor along the gene axis
    sliced_tensor = tensor_3d[selected_indices, :, :]

    return sliced_tensor, selected_genes


def cal_eigen_varexp(tensor, mode0, modes, varexp):
    """
    Perform eigen decomposition of the unfolded tensor's covariance matrix,
    plot the variance explained by each component, and return the number of
    components whose individual explained variance is greater than the threshold.

    Parameters:
    - tensor: a tensorly tensor object
    - mode0: the unfolding mode to use as rows
    - modes: the modes to use as columns (unused, kept for compatibility)
    - varexp: explained variance threshold (%) for individual components

    Returns:
    - fig: matplotlib figure object of the barplot
    - rank: number of components with explained variance > `varexp`%
    """
    # Unfold the tensor along the specified mode
    arr = unfold(tensor, mode=mode0)

    # Compute covariance-like matrix
    mat = np.dot(arr, arr.T)

    # Eigen decomposition (symmetric matrix, descending order)
    eigvals = np.linalg.eigvalsh(mat)[::-1]
    eig_prop = eigvals / np.sum(eigvals)  # normalized explained variance

    # Plot top 20 components
    top_k = min(20, len(eig_prop))
    ks = np.arange(1, top_k + 1)
    values = eig_prop[:top_k] * 100

    #fig, ax = plt.subplots(figsize=(10, 5))
    #ax.bar(ks, values, alpha=0.4)
    #ax.axhline(y=varexp, linestyle='dashed', color='red', label=f'Threshold {varexp}%')
    #ax.set_xticks(ks)
    #ax.set_xlabel("Component")
    #ax.set_ylabel("Explained Variance (%)")
    #ax.set_yscale("log")
    #ax.set_title("Explained Variance by Principal Components")
    #ax.legend()
    #plt.tight_layout()

    # Count how many components have explained variance > threshold
    threshold = varexp / 100.0
    rank = len(eig_prop) - np.searchsorted(eig_prop[::-1], threshold, side='right')

    return rank





def TensorReg_3d(tensor, X_covar1=None,core_shape=None):
    """
    Perform Tucker decomposition with covariates on a 3D tensor.

    Parameters:
    - tensor: numpy array or tensorly tensor of shape (d1, d2, d3)
    - X_covar1: sample-level covariates, shape (d1, p1); defaults to identity
    - core_shape: tuple (r1, r2, r3) for Tucker rank

    Returns:
    - dict with factor matrices W1–W3, core tensor G, reconstructed tensor U,
      intermediate C_ts, log-likelihood array lglk, and estimated variance sigma
    """
    d1, d2, d3 = tensor.shape
    r1, r2, r3 = core_shape

    # Use identity if covariate not provided
    if X_covar1 is None:
        X_covar1 = np.eye(d1)
    X_covar2 = np.eye(d2)
    X_covar3 = np.eye(d3)

    # QR decomposition of covariates
    Q1, R1 = np.linalg.qr(X_covar1)
    Q2, R2 = np.linalg.qr(X_covar2)
    Q3, R3 = np.linalg.qr(X_covar3)

    # Project tensor into orthogonal covariate space
    tensor_proj = multi_mode_dot(tensor, [Q1.T, Q2.T, Q3.T], modes=[0, 1, 2])


    # Project back to covariate space using R⁻¹
    R1_inv = np.linalg.pinv(R1)
    R2_inv = np.linalg.pinv(R2)
    R3_inv = np.linalg.pinv(R3)
    core, factors = tucker(tensor_proj, rank=core_shape)
    B_hat = tl.tucker_to_tensor((core, factors))

    C_ts = multi_mode_dot(B_hat, [R1_inv, R2_inv, R3_inv], modes=[0, 1, 2])

    # Orthogonal Tucker decomposition of C_ts
    core_final, factors = tucker(C_ts, rank=core_shape)
    W1, W2, W3 = factors
    G = core_final

    # Reconstruct the full tensor from core and covariates
    U = multi_mode_dot(C_ts, [X_covar1, X_covar2, X_covar3], modes=[0, 1, 2])

    # Residual variance and log-likelihood
    residual = tl.to_numpy(tensor) - U
    sigma_est = np.mean(residual**2)
    lglk = norm.logpdf(tl.to_numpy(tensor), loc=U, scale=np.sqrt(sigma_est))

    return {
        'W': {'W1': W1, 'W2': W2, 'W3': W3},
        'X_covar1': X_covar1,
        'G': G,
        'U': U,
        'C_ts': C_ts,
        'lglk': lglk,
        'sigma': sigma_est
    }



if __name__ == "__main__":

    # Dimensions
    n_samples, n_features2, n_features3 = 30, 20, 15
    Q = 10  # sample-level covariate dimension
    rq,r2, r3 = 10,4, 4  # Tucker rank

    # Sample-level covariates
    X_covar1 = np.random.randn(n_samples, Q)

    # Latent factor matrices (true W)
    W1 = np.random.randn(Q,rq)
    W2 = np.random.randn(n_features2, r2)
    W3 = np.random.randn(n_features3, r3)

    # Core tensor
    G = np.random.randn(rq, r2, r3)

    M_list = [W1,W2,W3]
    B = multi_mode_dot(G,M_list,modes=[0,1,2])
    Y = mode_dot(B, X_covar1, mode=0)
    noise = np.random.normal(scale=0.01, size=Y.shape)
    T_noisy = Y + noise
    
    res = TensorReg_3d(tensor=T_noisy,
                  X_covar1=X_covar1,
                  core_shape=(rq, r2, r3))