import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium", css_file="")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    return mo, np, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # שאלה 1

    1. תמונה להוסיף
    """)
    return


@app.cell
def _(mo):
    mo.image(src="t2q1.1.png",alt="להוסיף את התמונה של הזמנים אחד ביחס לשני")
    return


@app.cell
def _(mo):
    mo.md(r"""
    $$K_{xz1} < K_{xz2}$$

    $$K_{yz1} > K_{yz2}$$

    בתחילת האקטיבציה רמת ה
    X
    גבוהה יותר מרמת ה
    Y
    ולכן
    Z1
    מופעל מרמת ה
    X
    ו
    Z2
    מופעל מרמת ה
    Y
    שעוברים את קבוע הדיסאסוציאציה שלהם.

    בעצירה הסדר הוא הפוך וצריך ששני הגנים יהיו מתחת לקבוע הדיסאסוציאציה כדי להפסיק את הייצור,
    X
    דועך לפני
    Y
    ולכן הקבועים הפוכים בין ה
    Z
    ים השונים
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # שאלה 2

    א.
    """)
    return


@app.cell
def _():
    # setup for all
    #steps
    time = 100
    dt = 0.001
    steps = int(time / dt)

    # disasociation constants
    kd_xy = 1.2
    kd_xz = 1.2
    kd_yz = 1.2

    # X constants
    a_x = 0.000001 #dismanteling

    # Y constants
    a_y = 0.01

    # Z constants
    a_z = 1

    return a_x, a_y, a_z, dt, kd_xy, kd_xz, kd_yz, steps


@app.cell
def _(mo):
    mo.md(r"""
    $$\frac{dY}{dt} = \frac{X}{k_{d} + X} - \alpha_{y} Y$$

    $$\frac{dZ}{dt} = \frac{Y}{k_{dyz} +Y} - \alpha_{z}Z$$
    """)
    return


@app.cell
def _(a_y, a_z, dt, kd_xy, kd_xz, kd_yz, np, steps):
    ax = np.ones(steps) * kd_xy * 10
    ay = np.zeros(steps)
    az = np.zeros(steps)
    asx = np.zeros(steps)
    asx[1000:] = 1

    for _step in range(1, steps):
        _dy = ((ax[_step - 1]) / (ax[_step - 1] + kd_xy) - a_y * ay[_step - 1]) * dt if asx[_step] > 0 else - a_y * ay[_step - 1] * dt

        _dz = - a_z * az[_step - 1] * dt
        if ax[_step] >= kd_xz:
            _dz += (ax[_step - 1]) / (kd_xz + ax[_step - 1]) * dt
        if ay[_step] >= kd_yz:
            _dz += (ay[_step -1]) / (kd_yz + ay[_step -1]) * dt
        

        ay[_step] = ay[_step - 1] + _dy
        az[_step] = az[_step - 1] + _dz
    az
    return asx, ax, ay, az


@app.cell
def _(asx, ax, ay, az, plt, steps):
    plt.plot(range(steps), ax)
    plt.plot(range(steps), asx)
    plt.plot(range(steps), ay)
    plt.plot(range(steps), az)
    return


@app.cell
def _(a_x, a_y, a_z, ay, dt, kd_xy, kd_xz, kd_yz, np, steps):
    bx = np.zeros(steps)
    by = np.zeros(steps)
    bz = np.zeros(steps)
    bsx = np.zeros(steps)

    start_sx = 1000
    bsx[start_sx:] = 1

    for _i in range(start_sx, steps):
        bx[_i] = 1 - np.exp(- a_x * _i)

    for _step in range(1, steps):
    
        _dy = ((bx[_step - 1]) / (bx[_step - 1] + kd_xy) - a_y * by[_step - 1]) * dt if bsx[_step] > 0 else - a_y * by[_step - 1] * dt

        _dz = - a_z * bz[_step - 1] * dt
        if bx[_step] >= kd_xz:
            _dz += (bx[_step - 1]) / (kd_xz + bx[_step - 1]) * dt
        if by[_step] >= kd_yz:
            _dz += (by[_step -1]) / (kd_yz + ay[_step -1]) * dt
        

        by[_step] = by[_step - 1] + _dy
        bz[_step] = bz[_step - 1] + _dz
    bx[999:1010]
    return (bx,)


@app.cell
def _(bx, plt, steps):
    plt.plot(range(steps), bx)
    # plt.plot(range(steps), bsx)
    # plt.plot(range(steps), by)
    # plt.plot(range(steps), bz)
    return


if __name__ == "__main__":
    app.run()
