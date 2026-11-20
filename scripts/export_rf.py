# export_rf.py
import joblib

def export_tree(tree, f):
    tree_ = tree.tree_
    for i in range(tree_.node_count):
        if tree_.children_left[i] != -1:
            # internal node: feature_id threshold left_child right_child
            f.write(f"{tree_.feature[i]} {tree_.threshold[i]} {tree_.children_left[i]} {tree_.children_right[i]}\n")
        else:
            # leaf node: -1 0 0 0 leaf_value
            f.write(f"-1 0.0 0 0 {tree_.value[i][0][0]}\n")

def export_model(model, filename):
    with open(filename, 'w') as f:
        f.write(f"{model.n_estimators}\n") # first line: number of trees
        for estimator in model.estimators_:
            export_tree(estimator, f)

# load your model and export
model = joblib.load("rf_error_model.pkl")
export_model(model, "rf_model.txt")
print("model exported to rf_model.txt")