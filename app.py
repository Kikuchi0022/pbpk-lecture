import streamlit as st
st.set_page_config(
    page_title="PBPKモデル",
    layout="wide"
)

# ファイル名例: pk_model_app.py
import streamlit as st
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import font_manager
import matplotlib as mpl
import japanize_matplotlib


st.header("🧬 肝臓コンパートメントを含むPBPKモデル")
st.write("講義用")
# ユーザー入力
Dose = st.slider("投与量 (mg)", 10.0, 1000.0, 100.0, step=10.0)
Vb = st.slider("血漿容積 Vb (L)", 1.0, 10.0, 5.0, step=0.1)
VH = st.slider("肝臓容積 VH (L)", 0.5, 5.0, 1.5, step=0.1)
QH = st.slider("肝血流量 QH (L/h)", 0.1, 5.0, 1.2, step=0.1)
Rb = st.slider("血液/組織比 Rb", 0.5, 2.0, 1.0, step=0.1)
KPH = st.slider("肝分配係数 Kp,H", 0.5, 5.0, 1.0, step=0.1)
CLR = st.slider("腎クリアランス CL_R (L/h)", 0.0, 5.0, 0.5, step=0.1)
CLHint = st.slider("肝クリアランス CL_H,int (L/h)", 0.0, 5.0, 1.0, step=0.1)
fu_p = st.slider("非結合率 fu_p", 0.01, 1.0, 0.1, step=0.01)
ka = st.slider("吸収速度定数 ka (/h)", 0.1, 5.0, 1.0, step=0.1)
FaFg = st.slider("吸収率×腸管通過率 Fa×Fg", 0.1, 1.0, 0.8, step=0.05)

time_end = 24  # シミュレーション時間（h）
n_points = 200  # 時間分割数

# time_end = st.sidebar.slider("シミュレーション時間 (h)", 1, 48, 24)
# n_points = st.sidebar.slider("時間分割数", 50, 500, 200)




# 初期条件
Xg0 = FaFg * Dose
Cb0 = 0.0
CH0 = 0.0
y0 = [Cb0, CH0, Xg0]
t = np.linspace(0, time_end, n_points)

# 微分方程式
def model(y, t):
    Cb, CH, Xg = y
    dCb_dt = (-QH * Cb + QH * CH * Rb / KPH - CLR * Cb) / Vb
    dCH_dt = (QH * Cb - QH * CH * Rb / KPH + ka * Xg - CLHint * CH / KPH * fu_p) / VH
    dXg_dt = -ka * Xg
    return [dCb_dt, dCH_dt, dXg_dt]

# 数値解
sol = odeint(model, y0, t)
Cb, CH, Xg = sol.T

# グラフ描画
fig, ax = plt.subplots()
ax.plot(t, Cb, label="血漿濃度 Cb (mg/L)", color="blue")
ax.plot(t, CH, label="肝濃度 CH (mg/L)", color="red")
#ax.plot(t, Xg, label="消化管内量 Xg (mg)", color="green")
ax.set_xlabel("時間 (h)")
ax.set_ylabel("濃度 / 量")
ax.set_title("PBPKモデルによる濃度推移")
ax.legend()
ax.grid(True)
st.pyplot(fig)
