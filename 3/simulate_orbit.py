"""Простое моделирование орбиты по TLE CelesTrak с помощью SGP4."""

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sgp4.api import Satrec, WGS72
from sgp4.conveniences import sat_epoch_datetime


NORAD_ID = 57172
SECONDS_IN_DAY = 86400
DAYS_TO_SIMULATE = 7
STEP = 60
MU = 398600.8  # гравитационный параметр Земли, км³/с²
EARTH_RADIUS = 6378.135  # радиус Земли, км
FOLDER = Path(__file__).parent
DATA_FOLDER = FOLDER / "data"
RESULTS_FOLDER = FOLDER / "results"


def load_satellite():
    """Читает сохранённые двухстрочные элементы орбиты."""
    tle_file = DATA_FOLDER / f"tle_{NORAD_ID}.txt"
    with tle_file.open(encoding="utf-8") as file:
        name = file.readline().strip()
        line1 = file.readline().rstrip("\r\n")
        line2 = file.readline().rstrip("\r\n")

    satellite = Satrec.twoline2rv(line1, line2, WGS72)
    return name, satellite


def get_tle_elements(satellite):
    """Средние элементы, прочитанные библиотекой из TLE."""
    mean_motion = satellite.no_kozai / 60  # рад/с
    semi_major_axis = (MU / mean_motion**2) ** (1 / 3)
    eccentricity = satellite.ecco
    inclination = math.degrees(satellite.inclo)
    raan = math.degrees(satellite.nodeo)
    argument_of_perigee = math.degrees(satellite.argpo)
    return semi_major_axis, eccentricity, inclination, raan, argument_of_perigee


def get_elements(position, velocity):
    """Осцулирующие элементы по положению и скорости спутника."""
    r = np.array(position)
    v = np.array(velocity)
    r_length = np.linalg.norm(r)
    h = np.cross(r, v)
    h_length = np.linalg.norm(h)
    e_vector = np.cross(v, h) / MU - r / r_length
    eccentricity = np.linalg.norm(e_vector)
    energy = np.dot(v, v) / 2 - MU / r_length
    semi_major_axis = -MU / (2 * energy)

    inclination = math.degrees(math.acos(h[2] / h_length))
    node = np.array([-h[1], h[0], 0])
    node_length = np.linalg.norm(node)
    raan = math.degrees(math.atan2(node[1], node[0]))
    cosine = np.dot(node, e_vector) / (node_length * eccentricity)
    sine = np.dot(np.cross(node, e_vector), h) / (node_length * eccentricity * h_length)
    argument_of_perigee = math.degrees(math.atan2(sine, cosine))

    return semi_major_axis, eccentricity, inclination, raan, argument_of_perigee


def main():
    name, satellite = load_satellite()
    epoch = sat_epoch_datetime(satellite)
    times = np.arange(0, DAYS_TO_SIMULATE * SECONDS_IN_DAY + STEP, STEP)

    positions = []
    velocities = []
    orbit_elements = []

    for seconds in times:
        jd = satellite.jdsatepoch + satellite.jdsatepochF + seconds / SECONDS_IN_DAY
        whole_day = math.floor(jd)
        _, position, velocity = satellite.sgp4(whole_day, jd - whole_day)
        positions.append(position)
        velocities.append(velocity)
        orbit_elements.append(get_elements(position, velocity))

    positions = np.array(positions)
    velocities = np.array(velocities)
    orbit_elements = np.array(orbit_elements)
    radius = np.linalg.norm(positions, axis=1)
    speed = np.linalg.norm(velocities, axis=1)
    hours = times / 3600

    RESULTS_FOLDER.mkdir(exist_ok=True)

    # Траектория спутника вокруг Земли.
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    u = np.linspace(0, 2 * np.pi, 40)
    w = np.linspace(0, np.pi, 20)
    earth_x = EARTH_RADIUS * np.outer(np.cos(u), np.sin(w))
    earth_y = EARTH_RADIUS * np.outer(np.sin(u), np.sin(w))
    earth_z = EARTH_RADIUS * np.outer(np.ones_like(u), np.cos(w))
    ax.plot_surface(earth_x, earth_y, earth_z, color="skyblue", alpha=0.35)
    ax.plot(*positions.T, color="firebrick", linewidth=0.6)
    ax.set(xlabel="x, км", ylabel="y, км", zlabel="z, км", title=name)
    ax.set_box_aspect((1, 1, 1))
    fig.tight_layout()
    fig.savefig(RESULTS_FOLDER / "trajectory_3d.png", dpi=160)
    plt.close(fig)

    # Расстояние от центра Земли и скорость.
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
    axes[0].plot(hours, radius)
    axes[0].set_ylabel("r, км")
    axes[1].plot(hours, speed)
    axes[1].set_ylabel("v, км/с")
    axes[1].set_xlabel("Время от эпохи, ч")
    for axis in axes:
        axis.grid(True)
    fig.tight_layout()
    fig.savefig(RESULTS_FOLDER / "radius_speed.png", dpi=160)
    plt.close(fig)

    # Изменение осцулирующих элементов орбиты.
    titles = ["a, км", "e", "i, град", "r$\Omega$, град", r"$\omega$, град"]
    fig, axes = plt.subplots(5, 1, sharex=True, figsize=(9, 10))
    for column, axis, title in zip(range(5), axes, titles):
        values = orbit_elements[:, column]
        if column == 3:
            values = np.degrees(np.unwrap(np.radians(values)))
        axis.plot(hours, values)
        axis.set_ylabel(title)
        axis.grid(True)
    axes[-1].set_xlabel("Время от эпохи, ч")
    fig.tight_layout()
    fig.savefig(RESULTS_FOLDER / "osculating_elements.png", dpi=160)
    plt.close(fig)

    a, e, inclination, raan, argument_of_perigee = get_tle_elements(satellite)
    print(f"Эпоха TLE: {epoch.isoformat()}")
    print("Средние элементы TLE:")
    print(f"a = {a:.3f} км, e = {e:.7f}, i = {inclination:.4f} град")
    print(f"Omega = {raan:.4f} град, omega = {argument_of_perigee:.4f} град")
    print("Начальные осциллирующие элементы:")
    for title, value in zip(titles, orbit_elements[0]):
        print(f"  {title}: {value:.6f}")
    print(f"Расстояние: {radius.min():.3f}–{radius.max():.3f} км")
    print(f"Скорость: {speed.min():.5f}–{speed.max():.5f} км/с")

    print(f"{name}, NORAD {NORAD_ID}: {len(times)} точек, интервал {DAYS_TO_SIMULATE} суток.")
    print(f"Результаты сохранены в {RESULTS_FOLDER}")


if __name__ == "__main__":
    main()
