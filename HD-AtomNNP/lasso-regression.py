from sklearn import linear_model
import numpy as np

x = np.array([[0.0,1.0,2.0]
             ,[1.0, 2.0, 3.0]
             ,[2.0, 3.0, 4.0]
             ,[3.0, 4.0, 5.0]
             ,[4.0, 5.0, 6.0]
             ,[4.0, 4.0, 4.0]], dtype=np.float)

y = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 4.0], dtype=np.float)

clf = linear_model.Lasso(alpha=0.1, max_iter=1000)

clf.fit(x, y)

P = clf.predict([[3.0,4.0,5.0],[0.5,1.5,2.5]])

print (P)

print(clf.coef_)
print(clf.intercept_)