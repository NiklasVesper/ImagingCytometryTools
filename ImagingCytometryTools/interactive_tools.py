import cv2

def interactive_area_selection(image_path, protein, file_paths=None, additional_images=None):

    print("Press the 'r' to reset, press 'a' to add an additional area and 'c' to continue.")

    areas = []
    working_on_it = True

    while working_on_it == True:

        if file_paths != None:

            for index, file_path in enumerate(file_paths):
                add_images = cv2.imread(file_path, cv2.IMREAD_UNCHANGED)
                cv2.imshow('Additional image: ' + additional_images[index], add_images)

        image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
        cv2.imshow('Gating image: ' + protein, image)
        key = cv2.waitKey(0)

        if key == ord('c'):
            working_on_it = False

        while key == ord('a'):
            cv2.destroyWindow('Gating image: ' + protein)

            image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
            image = cv2.cvtColor(image, cv2.COLOR_GRAY2RGB)
            clone = image.copy()
            area = []

            def click_event(event, x, y, flags, param):
                if event == cv2.EVENT_LBUTTONDOWN:
                    area.append((x, y))
                    cv2.circle(image, (x, y), 3, (0, 0, 255), -1)
                    if len(area) > 1:
                        cv2.line(image, area[-2], area[-1], (255, 0, 0), 2)
                    cv2.imshow('Gating image: ' + protein, image)

            cv2.imshow('Gating image: ' + protein, image)
            cv2.setMouseCallback('Gating image: ' + protein, click_event)

            key = cv2.waitKey(0)

            if key == ord('r'):
                print('Resetted')
                cv2.destroyWindow('Gating image: ' + protein)
                continue

            if not area:
                continue
            else:
                areas.append(area)

            if key == ord('c'):
                working_on_it = False
                cv2.destroyWindow('Gating image: ' + protein)
                continue

        cv2.destroyAllWindows()

    cv2.destroyAllWindows()

    return(areas)

def interactive_compensation(x):
    return x