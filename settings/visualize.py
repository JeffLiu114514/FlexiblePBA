import sys
import json
import matplotlib.pyplot as plt

def vis(filename):
    settings = json.load(open(filename))
    plt.axvline(0, 0, 1, color='r')
    levels = []
    for level in settings["level_settings"]:
        plt.axvline(level["adc_upper_read_ref_lvl"], 0, 1, color='r')
        plt.axvline(level["adc_lower_write_ref_lvl"], 0, 0.5, color='g', linestyle='--')
        plt.axvline(level["adc_upper_write_ref_lvl"], 0, 0.5, color='g', linestyle='--')
        levels.append([level["adc_upper_read_ref_lvl"], level["adc_lower_write_ref_lvl"], level["adc_upper_write_ref_lvl"]])
    print("Upper read, lower write, upper write")
    print(levels)
    plt.xlim(-0.5, 63.5)
    plt.gca().get_yaxis().set_visible(False)
    plt.show()

if __name__ == "__main__":
    filename = sys.argv[1]
    vis(filename)
