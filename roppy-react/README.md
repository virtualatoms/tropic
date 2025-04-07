# Roppy React

A modern React-based web interface for the Roppy polymer analysis platform.

## Features

- Browse and search through a database of monomers
- View detailed information about each monomer
- Interactive 3D molecular visualization
- Responsive design with dark mode support
- RESTful API documentation

## Prerequisites

- Node.js (v16 or later)
- npm (v7 or later)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/roppy-react.git
cd roppy-react
```

2. Install dependencies:
```bash
npm install
```

3. Start the development server:
```bash
npm run dev
```

The application will be available at `http://localhost:5173`.

## Building for Production

To create a production build:

```bash
npm run build
```

The built files will be in the `dist` directory.

## Project Structure

```
src/
  ├── components/     # Reusable UI components
  ├── layouts/        # Layout components
  ├── pages/          # Page components
  ├── hooks/          # Custom React hooks
  ├── services/       # API service functions
  ├── types/          # TypeScript type definitions
  └── assets/         # Static assets
```

## Technologies Used

- React
- TypeScript
- Mantine UI
- 3Dmol.js
- React Router
- Vite

## License

MIT
