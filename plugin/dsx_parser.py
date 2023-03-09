import itertools


__all__ = ['dsx_parse']


def dsx_parse(dsx_output_path, ligand_comp):
    with open(dsx_output_path, 'r') as f:
        dsx_output_lines = f.readlines()

    scores = dict()
    ligand_lines = []
    while find_next_ligand(dsx_output_lines):
        ligand_lines = parse_ligand_lines(dsx_output_lines)
        for line in ligand_lines:
            line_items = line.split("__")
            atom_items = line_items[1].split("_")
            score = float(line_items[2])
            ireceptor_iligand = (int(atom_items[1]), int(atom_items[2]))
            if ireceptor_iligand not in scores:
                scores[ireceptor_iligand] = [score]
            else:
                scores[ireceptor_iligand].append(score)
        # Attach calculated score to each atom
        ligand_molecule = next(
            mol for i, mol in enumerate(ligand_comp.molecules)
            if i == ligand_comp.current_frame
        )
        if not hasattr(ligand_molecule, "atom_score_limits"):
            ligand_molecule.atom_score_limits = [float('inf'), float('-inf')]

        for atom_tuple, score_arr in scores.items():
            score = sum(score_arr) / len(score_arr)
            ligand_molecule.atom_score_limits[0] = min(ligand_molecule.atom_score_limits[0], score)
            ligand_molecule.atom_score_limits[1] = max(ligand_molecule.atom_score_limits[1], score)
            atom = next(itertools.islice(ligand_molecule.atoms, atom_tuple[1] - 1, atom_tuple[1]))
            atom.score = score


def parse_ligand_lines(dsx_output_lines):
    ligand_lines = []
    endline = "# End of pair potentials"
    while True:
        line = dsx_output_lines.pop(0)
        ligand_lines.append(dsx_output_lines.pop(0))
        if line.startswith(endline):
            # Pop empty line with '----------\n'
            ligand_lines.pop(-1)
            break
    return ligand_lines


def find_next_ligand(line_list):
    ligand_line_start = "# Receptor-Ligand:"
    starting_line_index = 0
    ligand_found = False
    for i in range(0, len(line_list)):
        line = line_list[i]
        if line.startswith(ligand_line_start):
            starting_line_index = i
            ligand_found = True
            break
    del line_list[:starting_line_index]
    return ligand_found


def parse_scores(self, dsx_output):
    lines = dsx_output.splitlines()
    number_of_lines = len(lines)
    line_index = 0
    if not find_next_ligand():
        self.stop_scoring()
        self.clear_sphere_streams()
        return

    ligand_index = 0
    has_next_ligand = True
    while has_next_ligand:
        scores = dict()
        last_tuple = None
        last_arr = None
        score_min = None
        score_max = None

        while line_index < number_of_lines - 1:
            line = lines[line_index]
            if line.startswith("# End of pair potentials"):
                has_next_ligand = find_next_ligand()
                break
            line_items = line.split("__")
            atom_items = line_items[1].split("_")
            score = float(line_items[2])
            ireceptor_iligand = (int(atom_items[1]), int(atom_items[2]))

            if last_tuple != ireceptor_iligand:
                if ireceptor_iligand in scores:
                    last_arr = scores[ireceptor_iligand]
                else:
                    last_arr = []
                    scores[ireceptor_iligand] = last_arr
            last_tuple = ireceptor_iligand
            last_arr.append(score)

            if score_min is None or score < score_min:
                score_min = score
            if score_max is None or score > score_max:
                score_max = score
            line_index += 1

        if score_max is None or score_min is None:
            continue
        # score_gap = max(score_max - score_min, 0.01)
        try:
            for atom_tuple, score_arr in scores.items():
                score = sum(score_arr) / len(score_arr)
                # bfactor = ((score - score_min) / score_gap) * BFACTOR_GAP + BFACTOR_MIN
                # bfactor_two = score / (-score_min if score < 0 else score_max)
                molecule = self._ligands._molecules[atom_tuple[0] - 1 + ligand_index]
                if not hasattr(molecule, "atom_score_limits"):
                    molecule.atom_score_limits = [float('inf'), float('-inf')]
                if score < molecule.atom_score_limits[0]:
                    molecule.atom_score_limits[0] = score
                elif score > molecule.atom_score_limits[1]:
                    molecule.atom_score_limits[1] = score

                atom = next(itertools.islice(molecule.atoms, atom_tuple[1] - 1, atom_tuple[1]))
                atom.score = score
        except Exception:
            err_msg = "Error parsing ligand scores. Are your ligands missing bonds?"
            # self.send_notification(NotificationTypes.error, err_msg)
            # nanome.util.Logs.error(err_msg)
            self.stop_scoring()
            return

        ligand_index += 1

    labels = []
    colors = []
    scales = []
    for atom in self._ligands.atoms:
        red = 0
        green = 0
        blue = 0
        alpha = 141
        if hasattr(atom, "score"):
            ligand = atom.molecule
            denominator = -ligand.atom_score_limits[0] if atom.score < 0 else ligand.atom_score_limits[1]
            norm_score = atom.score / denominator
            red = 255 if norm_score > 0 else 0
            green = 255 if -norm_score >= 0 else 0
            red_scale = int(max(norm_score * 255, 0))
            green_scale = int(max(-norm_score * 255, 0))
            scale = max(green_scale, red_scale) / 255. + 1
        else:
            norm_score = 0.0
            green = 0
            red = 255
            scale = 1

        labels.append(f'{norm_score:.3f}')
        colors.append(red)
        colors.append(green)
        colors.append(blue)
        colors.append(alpha)
        scales.append(scale)

    try:
        if self._streams_ready:
            if self._label_stream:
                self._label_stream.update(labels)
            self._color_stream.update(colors)
            self._scale_stream.update(scales)
    except Exception:
        print("Error while updating sphere stream")
