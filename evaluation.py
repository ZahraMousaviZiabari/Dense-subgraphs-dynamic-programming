import numpy as np
from scipy import stats

def perform_hypothesis_test_rel(real_data, synthetic_data):
    # Assume real_data and synthetic_data are lists of performance metrics
    # obtained from the real dynamic network and synthetic graph, respectively.

    # Perform a two-sided t-test
    t_statistic, p_value = stats.ttest_rel(real_data, synthetic_data)

    return t_statistic, p_value

# # Example data (replace with actual algorithm performance metrics)
# real_data = [0.85, 0.92, 0.88, 0.91, 0.89]
# synthetic_data = [0.78, 0.80, 0.79, 0.82, 0.81]




def perform_hypothesis_test_ind(real_data, synthetic_data):
    # Assume real_data and synthetic_data are lists of lists
    # where each inner list represents the performance metrics for a condition.

    # Perform a two-sided independent-sample t-test
    t_statistic, p_value = stats.ttest_ind(real_data, synthetic_data)

    return t_statistic, p_value

# # Example data (replace with actual algorithm performance metrics)
# real_data = [[0.85, 0.92, 0.88, 0.91, 0.89],
#              [0.90, 0.88, 0.87, 0.91, 0.92]]

# synthetic_data = [[0.78, 0.80, 0.79, 0.82, 0.81],
#                   [0.75, 0.79, 0.81, 0.78, 0.80, 0.85]]  # Different number of observations

    # print("T-statistic:", t_statistic)
    # print("P-value:", p_value)

    # alpha = 0.05  # Significance level
    # if p_value < alpha:
    #     print("Reject the null hypothesis. There is a significant difference.")
    # else:
    #     print("Fail to reject the null hypothesis. No significant difference.")


    
def calculate_z_scores(real_data, synthetic_data):
    
    # Calculate mean and standard deviation
    mean_real = np.mean(real_data) 
    std_dev_real = np.std(real_data)
    
    # Calculate z-scores
    z_scores_real = (real_data - mean_real) / std_dev_real
    z_scores_synthetic = (synthetic_data - mean_real) / std_dev_real
    
    avg_z_scores_real = np.mean(z_scores_real)
    avg_z_scores_synthetic = np.mean(z_scores_synthetic)
    
    print("Z-scores for real data:", z_scores_real) 
    print("Average Z-scores for real data:", avg_z_scores_real) 
    print("Z-scores for synthetic data:", z_scores_synthetic) 
    print("Average Z-scores for synthetic data:", avg_z_scores_synthetic) 
    return

