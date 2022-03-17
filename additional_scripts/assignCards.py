import json
import shutil
import os

with open("card_assignment.json") as f:
 card_assignment = json.load(f)

for proc, cards in card_assignment.items():
  shutil.copyfile(os.path.join("cards", "run_cards", cards[0]), os.path.join("cards", proc, "run_card.dat"))
  shutil.copyfile(os.path.join("cards", "pythia_cards", cards[1]), os.path.join("cards", proc, "pythia8_card.dat"))
