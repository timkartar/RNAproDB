import requests
import os


def download_assembly_images(pdb_ids, output_folder="../rnaprodb_frontend/public/pdb_thumbnails"):
    """
    Downloads Assembly 1 images for a list of PDB IDs.
    """

    for pdb_id in pdb_ids:
        image_url = f"https://cdn.rcsb.org/images/structures/{pdb_id}_assembly-1.jpeg"
        response = requests.get(image_url)
        if response.status_code == 200:
            with open(os.path.join(output_folder, f"{pdb_id}_assembly1.png"), "wb") as file:
                file.write(response.content)
            os.chmod(os.path.join(output_folder, f"{pdb_id}_assembly1.png"), 777)
        else:
            print(f"Failed to download image for {pdb_id}")

def main():
    with open("nakb_prna_ids.txt", "r") as file:
        pdb_ids = file.read().split(",")

    download_assembly_images(pdb_ids)

if __name__ == "__main__":
    main()
