
# EDQNM MHD alpha D3 MODEL 

## Description of the model
 - This code solves the 3-dimensional Eddy Damped Quasi-Normal Markovian Magneto Hydrodynamical (EDQNM-MHD) equation, refer to *Orsag 1970, Pouquet 1976* for more details.
 - In this model, a variation is implemented in the eddy damping timescale that is parameterised by $\alpha$, hence the name.

## Numerical simulation setup
- This is a one-dimensional system, described in 'N' wave-numbers, which are equally spaced in logarithmic scale.
$$ k_i = k_0 \lambda ^{(i-1)} \;, \; \{i=1,2,\cdots,N\} $$
Every wavenumber $k_i$ is associated with a band of width $\Delta k_i$ from $k_i^{(-)}$ to $k_i^{(+)}$, given by:
$$\Delta k_i=k_i \ln(\lambda) \;,\; k_i^{(-)}=\frac{\Delta k_i}{\lambda -1}\; , \; k_i^{(+)}=\frac{\lambda \Delta k_i}{\lambda-1} $$
In this discrete setting, the calculations involving integration is approximated by:
$$ \int_0^{\infty} dk \: g(k) \approxeq \sum_{i=1}^{N} g(k_i)\Delta k_i$$
- Generally not *recommended* to change the spacing $\lambda = 2^{1/4}$ and the base wave-number $k_0=2^{-3}$ from their set values, unless necessary. With these parameters, $N=53$ yields a $k_{\mathsf{max}}=k_N=1024$ and doubles when $N$ increases by $4$.

### Inputs for the simulation 
1. Edit **input.dat** to choose the size of the system $N$, evolution time $T$, and the number of times the data has to be saved $S$ in uniform time intervals.

> N - Shells
37
T - Total time of the simulation
10.0
S - No of saves to be made
20

2. Edit **system_basicvariables.f90** to further tune the simulation further:

		visc_status                            = 1
		! '1' TO INCLUDE VISCOSITY, '0' FOR INVISCID CASE
		
		diff_status                            = 1
		! '1' TO INCLUDE DIFFUSIVITY, '0' FOR IDEAL CASE

		forc_status                            = 1
		! '1' TO ACTIVATE FORCING, '0' FOR DECAYING CASE (ONLY FOR KINETIC SPECTRUM)

		coupling_status                        = 1
		! '1' FOR COUPLED CASE
		! '0' FOR PURE KINETIC CASE
		! '2' FOR ONE-WAY COUPLED CASE, WITH A FIXED E(K)

3. Choose $\alpha$ here:

		eddy_damping_exp                       = two/thr
		! THE EXPONENT IN THE MODEL

4. Viscosity $\nu$ and Diffusivity $\zeta$ are chosen automatically from the chosen resolution $N$, as the minimum values the system can model without thermalization. For custom values, edit the file as follows:

		! visc                                 = 0.01D0
		! UNCOMMENT FOR CUSTOM VISCOSITY

		! diff                                 = 0.02D0
		! UNCOMMENT FOR CUSTOM DIFFUSIVITY 

5. The total energy $\mathcal{E}$, kinetic $\mathcal{E}^{(u)}$ and magnetic $\mathcal{E}^{(b)}$ energy at $t=0$ is set here:


		energy_0                               = one
		! TOTAL ENERGY 

		energy_B_0                             = 1E-2
		! INITIAL MAGNETIC ENERGY

		energy_V_0                             = energy_0 - energy_B_0
		! INITIAL KINETIC ENERGY

6. The time-step is chosen with the aid of a CFL number and a smallest time-scale in the system. The rms time-scales and the viscous (diffusive) timescales are, respectively: 
$$ t_{\mathsf{rms}}^{(b)} = 1/\sqrt{\mathcal{E}^{(b)} k_{\mathsf{max}}^2} \quad, \quad t_{\mathsf{rms}}^{(u)} = 1/\sqrt{\mathcal{E}^{(u)} k_{\mathsf{max}}^2}$$
$$ \tau^{(b)} = ( \zeta k_{\mathsf{max}}^2 )^{-1} \quad, \quad \tau^{(u)} = ( \nu k_{\mathsf{max}}^2 )^{-1}  $$
Now for a given CFL (must be $>1$), the time step $\Delta t$ is chosen as:

		cfl_ref                                = 20
		! RATIO BETWEEN SMALLEST TIME SCALE IN THE SYSTEM AND THE TIME-STEP
$$ \Delta t = \frac{ \mathrm{min}: \{t^{(b)}_{\mathsf{rms}},t^{(u)}_{\mathsf{rms}}, \tau^{(b)},\tau^{(u)} \}}{\mathsf{CFL}}$$

### Forcing for the kinetic spectrum 
The forcing $F(k)$ concentrated mainly at the small wave-numbers is of the form:
$$ F(k_i)= A \; \epsilon \; k_i^2 \exp\left( -\frac{ k_i^2 } { 2k_I^2 }  \right) \; , \; \sum_{i=0}^{N}F(k_i)\Delta k_i = \epsilon $$	

### Small-scale dynamo perturbation 
The perturbation in the magnetic energy spectrum $E^{(b)}(k)$ at $t=0$, to investigate the small-scale dynamo scales as:
$$
E^{(b)}(k_i,0) \sim \exp \left( -\frac{ (k_i-k_{D})^2 }{k_0^2} \right)
$$


# Welcome to StackEdit!

Hi! I'm your first Markdown file in **StackEdit**. If you want to learn about StackEdit, you can read me. If you want to play with Markdown, you can edit me. Once you have finished with me, you can create new files by opening the **file explorer** on the left corner of the navigation bar.


# Files

StackEdit stores your files in your browser, which means all your files are automatically saved locally and are accessible **offline!**

## Create files and folders

The file explorer is accessible using the button in left corner of the navigation bar. You can create a new file by clicking the **New file** button in the file explorer. You can also create folders by clicking the **New folder** button.

## Switch to another file

All your files and folders are presented as a tree in the file explorer. You can switch from one to another by clicking a file in the tree.

## Rename a file

You can rename the current file by clicking the file name in the navigation bar or by clicking the **Rename** button in the file explorer.

## Delete a file

You can delete the current file by clicking the **Remove** button in the file explorer. The file will be moved into the **Trash** folder and automatically deleted after 7 days of inactivity.

## Export a file

You can export the current file by clicking **Export to disk** in the menu. You can choose to export the file as plain Markdown, as HTML using a Handlebars template or as a PDF.


# Synchronization

Synchronization is one of the biggest features of StackEdit. It enables you to synchronize any file in your workspace with other files stored in your **Google Drive**, your **Dropbox** and your **GitHub** accounts. This allows you to keep writing on other devices, collaborate with people you share the file with, integrate easily into your workflow... The synchronization mechanism takes place every minute in the background, downloading, merging, and uploading file modifications.

There are two types of synchronization and they can complement each other:

- The workspace synchronization will sync all your files, folders and settings automatically. This will allow you to fetch your workspace on any other device.
	> To start syncing your workspace, just sign in with Google in the menu.

- The file synchronization will keep one file of the workspace synced with one or multiple files in **Google Drive**, **Dropbox** or **GitHub**.
	> Before starting to sync files, you must link an account in the **Synchronize** sub-menu.

## Open a file

You can open a file from **Google Drive**, **Dropbox** or **GitHub** by opening the **Synchronize** sub-menu and clicking **Open from**. Once opened in the workspace, any modification in the file will be automatically synced.

## Save a file

You can save any file of the workspace to **Google Drive**, **Dropbox** or **GitHub** by opening the **Synchronize** sub-menu and clicking **Save on**. Even if a file in the workspace is already synced, you can save it to another location. StackEdit can sync one file with multiple locations and accounts.

## Synchronize a file

Once your file is linked to a synchronized location, StackEdit will periodically synchronize it by downloading/uploading any modification. A merge will be performed if necessary and conflicts will be resolved.

If you just have modified your file and you want to force syncing, click the **Synchronize now** button in the navigation bar.

> **Note:** The **Synchronize now** button is disabled if you have no file to synchronize.

## Manage file synchronization

Since one file can be synced with multiple locations, you can list and manage synchronized locations by clicking **File synchronization** in the **Synchronize** sub-menu. This allows you to list and remove synchronized locations that are linked to your file.


# Publication

Publishing in StackEdit makes it simple for you to publish online your files. Once you're happy with a file, you can publish it to different hosting platforms like **Blogger**, **Dropbox**, **Gist**, **GitHub**, **Google Drive**, **WordPress** and **Zendesk**. With [Handlebars templates](http://handlebarsjs.com/), you have full control over what you export.

> Before starting to publish, you must link an account in the **Publish** sub-menu.

## Publish a File

You can publish your file by opening the **Publish** sub-menu and by clicking **Publish to**. For some locations, you can choose between the following formats:

- Markdown: publish the Markdown text on a website that can interpret it (**GitHub** for instance),
- HTML: publish the file converted to HTML via a Handlebars template (on a blog for example).

## Update a publication

After publishing, StackEdit keeps your file linked to that publication which makes it easy for you to re-publish it. Once you have modified your file and you want to update your publication, click on the **Publish now** button in the navigation bar.

> **Note:** The **Publish now** button is disabled if your file has not been published yet.

## Manage file publication

Since one file can be published to multiple locations, you can list and manage publish locations by clicking **File publication** in the **Publish** sub-menu. This allows you to list and remove publication locations that are linked to your file.


# Markdown extensions

StackEdit extends the standard Markdown syntax by adding extra **Markdown extensions**, providing you with some nice features.

> **ProTip:** You can disable any **Markdown extension** in the **File properties** dialog.


## SmartyPants

SmartyPants converts ASCII punctuation characters into "smart" typographic punctuation HTML entities. For example:

|                |ASCII                          |HTML                         |
|----------------|-------------------------------|-----------------------------|
|Single backticks|`'Isn't this fun?'`            |'Isn't this fun?'            |
|Quotes          |`"Isn't this fun?"`            |"Isn't this fun?"            |
|Dashes          |`-- is en-dash, --- is em-dash`|-- is en-dash, --- is em-dash|


## KaTeX

You can render LaTeX mathematical expressions using [KaTeX](https://khan.github.io/KaTeX/):

The *Gamma function* satisfying $\Gamma(n) = (n-1)!\quad\forall n\in\mathbb N$ is via the Euler integral

$$
\Gamma(z) = \int_0^\infty t^{z-1}e^{-t}dt\,.
$$

> You can find more information about **LaTeX** mathematical expressions [here](http://meta.math.stackexchange.com/questions/5020/mathjax-basic-tutorial-and-quick-reference).


## UML diagrams

You can render UML diagrams using [Mermaid](https://mermaidjs.github.io/). For example, this will produce a sequence diagram:

```mermaid
sequenceDiagram
Alice ->> Bob: Hello Bob, how are you?
Bob-->>John: How about you John?
Bob--x Alice: I am good thanks!
Bob-x John: I am good thanks!
Note right of John: Bob thinks a long<br/>long time, so long<br/>that the text does<br/>not fit on a row.

Bob-->Alice: Checking with John...
Alice->John: Yes... John, how are you?
```

And this will produce a flow chart:

```mermaid
graph LR
A[Square Rect] -- Link text --> B((Circle))
A --> C(Round Rect)
B --> D{Rhombus}
C --> D
```